#' @title rciplot
#'
#' @description Create a scatterplot of your sample in which the x-axis maps to
#' the pre-scores, the y-axis maps to the post-scores and several graphical
#' elements (lines, colors) allow you to gain a quick overview about reliable
#' changes in these scores.
#' An example of this kind of plot is Figure 2 of Jacobson & Truax (1991).
#' Jacobson-Truax classification (represented in point colors) is always based
#' on `recovery_cutoff`, not on any other plotted horizontal line (e.g. mid of
#' means).
#'
#' @param data Dataframe containing all relevant data
#' @param pre Name of the column in `data` containing pre values
#' @param post Name of the column in `data` containing post values
#' @param group Name of column by which cases are to be grouped (controls shape
#' of scatter plot points)
#' @param reliability Reliability of the used test / instrument
#' @param recovery_cutoff Test score below which individuals are considered
#' healthy / recovered
#' @param show_classification_counts If TRUE, show number of cases for each
#' classification (e.g. reliable improvement, no reliable change, ...) in legend
#' @param show_classification_percentages Expanding on
#' `show_classification_counts`.If TRUE, show the respective percentage of the
#' whole sample each classification makes up.
#' @param higher_is_better TRUE if higher values indicate a remission / healthy
#' individual. FALSE if higher values indicate worse health.
#' @param pre_jitter Jitter factor to apply to pre values
#' @param post_jitter Jitter factor to apply to post values
#' @param opacity Alpha value of scatter plot points
#' @param size_points Size of scatter plot points.
#' @param size_lines Size (thickness) of lines in plot.
#' @param draw_meanmid_line Draw a horizontal line indicating the middle between
#' the population means for a functional (healthy) population and a
#' dysfunctional (diseased) population, described as criterion *c* in Jacobson &
#' Truax (1991).
#' @param draw_2sd_functional_line Draw a horizontal line indicating a cutoff
#' at a 2 SD distance from `mean_functional`, described as criterion *b* in
#' Jacobson & Truax (1991).
#' @param draw_2sd_dysfunctional_line Draw a horizontal line indicating a cutoff
#' at a 2 SD distance from `mean_dysfunctional`, described as criterion *a* in
#' Jacobson & Truax (1991).
#' @param mean_functional Required if `draw_meanmid_line = T` or
#' `draw_2sd_[dys]functional_line = T`.
#' Mean test score of the functional population.
#' @param mean_dysfunctional Required if `draw_meanmid_line = T` or
#' `draw_2sd_[dys]functional_line`.
#' Mean test score of the dysfunctional population.
#' @param sd_functional Optional for `draw_meanmid_line = T`. Standard deviation
#' of the functional population.
#' @param sd_dysfunctional Optional for `draw_meanmid_line = T`. Standard
#' deviation of the dysfunctional population.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{higher_is_better} \tab
#'        Exactly the input parameter \code{higher_is_better}
#'        \cr
#'    \tab \cr
#'    \code{reliable_change} \tab
#'        Pre-Post differences larger than this difference are regarded reliable
#'        \cr
#'    \tab \cr
#'    \code{plot} \tab
#'        ggplot2 scatter plot analogous to Figure 2 of Jacobson & Truax (1991)
#'        \cr
#' }
#'
#' @examples
#' # Using example data from `sample_data.rda` to recreate Figure 2 of
#' # Jacobson & Truax (1991):
#' rciplot(
#'     data = sample_data,
#'     pre = 'pre_data',
#'     post = 'post_data',
#'     reliability = 0.88,
#'     recovery_cutoff = 104,
#'     opacity = 1
#' )
#'
#' @importFrom dplyr case_when filter mutate mutate_at pull recode
#' @import ggplot2
#' @importFrom stats sd
#' @import tibble
#' @export

rciplot <- function(
    data,
    pre = NULL,
    post = NULL,
    group = NULL,
    reliability = NULL,
    recovery_cutoff = NULL,
    show_classification_counts = TRUE,
    show_classification_percentages = TRUE,
    higher_is_better = TRUE,
    pre_jitter = 0,
    post_jitter = 0,
    opacity = 0.5,
    size_points = 1,
    size_lines = 0.3,
    draw_meanmid_line = FALSE,
    draw_2sd_functional_line = FALSE,
    draw_2sd_dysfunctional_line = FALSE,
    mean_functional = NULL,
    mean_dysfunctional = NULL,
    sd_functional = 1,
    sd_dysfunctional = 1
) {
    if (
        is.null(pre) |
        is.null(post) |
        is.null(reliability)
    ) {
        stop(
            'These parameters are not optional: pre, post, reliability',
            call. = F
        )
    }

    pre_c <- data %>%
        pull(pre)
    post_c <- data %>%
        pull(post)

    if (pre_jitter)
        pre_c <- jitter(pre_c, pre_jitter)
    if (post_jitter)
        post_c <- jitter(post_c, post_jitter)

    # Calc Critical Difference ---------------------------------------
    # RC = (x2 -x1) / S_diff
    # S_diff = sqrt(2 * (SEM)^2)
    # SEM = sd * sqrt(1 - rel)
    # ... with:
    #   RC        = Reliable Change (pre-post-difference adjusted for
    #               reliability of measurement)
    #   S_diff    = Standard Error of the Difference (... of the
    #               pre-post-difference, of course)
    #   SEM       = Standard Error of Measurement (as estimated based on given
    #               standard deviation and reliability)
    #   sd        = Standard Deviation of the test
    #   rel       = Reliability (ideally retest-reliability) of the test
    # ... in which i want to solve the first equation for RC = 1.96 to find
    #   the minimum reliable difference x1-x2 which is going to provide the
    #   cutoff bands in the graph:
    #   1.96 = (x2 - x1) / S_diff  | * S_diff
    #   1.96 * S_diff = (x2 - x1)  | solved!
    sem = sd(pre_c) * sqrt(1 - reliability)
    s_diff = sqrt(2 * sem^2)
    diff_crit = 1.96 * s_diff

    # Draw plot ------------------------------------------------------

    # Build `df_plot` ......................................
    df_plot <- tibble(
        col_pre  = pre_c,
        col_post = post_c,
    )
    if (!is.null(group))
        df_plot$group <- data %>% pull(group)
    df_plot$diff <- df_plot$col_post - df_plot$col_pre
    if (higher_is_better) {
        df_plot <- df_plot %>%
            mutate(category = case_when(
                (diff > diff_crit) & col_post > recovery_cutoff ~
                    'reliable recovery',
                (diff > diff_crit) ~
                    'reliable improvement',
                (diff < (diff_crit * -1)) ~
                    'reliable deterioration',
                (col_post > recovery_cutoff) ~
                    'non-reliable recovery',
                T ~
                    'no reliable change'
            )) %>%
            mutate_at('category', as.factor)
    } else {
        df_plot <- df_plot %>%
            mutate(category = case_when(
                (diff < (diff_crit * -1)) & col_post < recovery_cutoff ~
                    'reliable recovery',
                (diff < (diff_crit * -1)) ~
                    'reliable improvement',
                (diff > diff_crit) ~
                    'reliable deterioration',
                (col_post < recovery_cutoff) ~
                    'non-reliable recovery',
                T ~
                    'no reliable change'
            )) %>%
            mutate_at('category', as.factor)
    }

    if (show_classification_counts) {
        count_relrec <- df_plot %>%
            dplyr::filter(category == 'reliable recovery') %>%
            nrow()
        count_relimp <- df_plot %>%
            dplyr::filter(category == 'reliable improvement') %>%
            nrow()
        count_reldet <- df_plot %>%
            dplyr::filter(category == 'reliable deterioration') %>%
            nrow()
        count_nonrec <- df_plot %>%
            dplyr::filter(category == 'non-reliable recovery') %>%
            nrow()
        count_nochan <- df_plot %>%
            dplyr::filter(category == 'no reliable change') %>%
            nrow()
        label_relrec <- paste0(
            'reliable recovery: ',
            count_relrec
        )
        label_relimp <- paste0(
            'reliable improvement: ',
            count_relimp
        )
        label_reldet <- paste0(
            'reliable deterioration: ',
            count_reldet
        )
        label_nonrec <- paste0(
            'non-reliable recovery: ',
            count_nonrec
        )
        label_nochan <- paste0(
            'no reliable change: ',
            count_nochan
        )
        if (show_classification_percentages) {
            label_relrec <- paste0(
                label_relrec,
                ' (', round(count_relrec / nrow(data), digits=4) * 100, '%)'
            )
            label_relimp <- paste0(
                label_relimp,
                ' (', round(count_relimp / nrow(data), digits=4) * 100, '%)'
            )
            label_reldet <- paste0(
                label_reldet,
                ' (', round(count_reldet / nrow(data), digits=4) * 100, '%)'
            )
            label_nonrec <- paste0(
                label_nonrec,
                ' (', round(count_nonrec / nrow(data), digits=4) * 100, '%)'
            )
            label_nochan <- paste0(
                label_nochan,
                ' (', round(count_nochan / nrow(data), digits=4) * 100, '%)'
            )
        }
        df_plot <- df_plot %>%
            mutate(category = recode(
                category,
                'reliable recovery' = label_relrec,
                'reliable improvement' = label_relimp,
                'reliable deterioration' = label_reldet,
                'non-reliable recovery' = label_nonrec,
                'no reliable change' = label_nochan,
            ))
    }

    # Draw plot ............................................
    totalmin = min(
        min(pre_c),
        min(post_c)
    )
    totalmax = max(
        max(pre_c),
        max(post_c)
    )
    # Scatter plot ///////////////////////////////
    plot <- df_plot %>%
        ggplot(aes(
            x = col_pre,
            y = col_post,
            color = category
        )) +
        geom_point(
            aes(
                shape = group
            ),
            alpha = opacity,
            size = size_points,
        ) +
        scale_color_manual(values=c(
            '#3891A6',  # no reliable change
            '#DEC402',  # non-reliable recovery
            '#DB5461',  # reliable deterioration
            '#CCD263',  # reliable improvement
            '#9BBC79'   # reliable recovery
        )) +
        scale_x_continuous(limits = c(totalmin, totalmax)) +
        scale_y_continuous(limits = c(totalmin, totalmax))

    if (!is.null(group)) {
        plot <- plot +
            scale_shape_manual(values=c(1, 4))
    }

    # Cutoff lines ///////////////////////////////
    df_lines = tibble(
        slope =     c(1,          1,         1,               0),
        intercept = c(0, -diff_crit, diff_crit, recovery_cutoff),
        name = c(
            'no change',
            'reliable change boundary',
            'reliable change boundary',
            paste0('recovery cutoff (', recovery_cutoff, ')')
        )
    )
    linetype_legend_order = c(
        'no change',
        'reliable change boundary',
        paste0('recovery cutoff (', recovery_cutoff, ')')
    )
    if (draw_meanmid_line) {
        midmean <- (
            (
                (sd_functional * mean_functional) +
                (sd_dysfunctional * mean_dysfunctional)
            ) / (
                sd_functional + sd_dysfunctional
            )
        )
        meanmid_line_name <- paste0('mid of means criterion (', midmean, ')')
        df_lines <- df_lines %>%
            add_row(
                slope = 0,
                intercept = midmean,
                name = meanmid_line_name
            )
        linetype_legend_order <- append(
            linetype_legend_order,
            meanmid_line_name
        )
    }
    if (draw_2sd_dysfunctional_line) {
        twosd_dysfunctional_cutoff <- mean_dysfunctional
        if (higher_is_better) {
            twosd_dysfunctional_cutoff <- twosd_dysfunctional_cutoff +
                (2 * sd_dysfunctional)
        } else {
            twosd_dysfunctional_cutoff <- twosd_dysfunctional_cutoff -
                (2 * sd_functional)
        }
        twosd_dysfunctional_name <- paste0(
            '2 SD from dysfunctional mean (',
            twosd_dysfunctional_cutoff,
            ')'
        )
        df_lines <- df_lines %>%
            add_row(
                slope = 0,
                intercept = twosd_dysfunctional_cutoff,
                name = twosd_dysfunctional_name
            )
        linetype_legend_order <- append(
            linetype_legend_order,
            twosd_dysfunctional_name
        )
    }
    if (draw_2sd_functional_line) {
        twosd_functional_cutoff <- mean_functional
        if (higher_is_better) {
            twosd_functional_cutoff <- twosd_functional_cutoff -
                (2 * mean(sd_functional, sd_dysfunctional))
        } else {
            twosd_functional_cutoff <- twosd_functional_cutoff +
                (2 * mean(sd_functional, sd_dysfunctional))
        }
        twosd_functional_name <- paste0(
            '2 SD from functional mean (',
            twosd_functional_cutoff,
            ')'
        )
        df_lines <- df_lines %>%
            add_row(
                slope = 0,
                intercept = twosd_functional_cutoff,
                name = twosd_functional_name
            )
        linetype_legend_order <- append(
            linetype_legend_order,
            twosd_functional_name
        )
    }
    plot <- plot +
        geom_abline(
            data = df_lines,
            mapping = aes(
                slope = slope,
                intercept = intercept,
                linetype = factor(name)
            ),
            alpha = 0.5,
            size = size_lines
        ) +
        scale_linetype_manual(
            values = c(
                'solid',
                'dashed',
                'dotted',
                'dotdash',
                'twodash',
                'longdash',
                'F1',
                '12345678'
            ),
            breaks = linetype_legend_order
        )

    # Theme & legend /////////////////////////////
    plot <- plot +
        theme_bw(base_size=14) +
        theme(
            legend.position='right',
        ) +
        labs(
            title = 'Jacobson-Truax plot',
            x = 'pre',
            y = 'post',
            color = 'Jacobson-Truax\nclassification',
            shape = 'Group',
            linetype = ''
        ) +
        guides(
            color = guide_legend(order = 1),
            shape = guide_legend(shape = 1),
            linetype = guide_legend(linetype = 1),
        )

    output <- list(
        higher_is_better = higher_is_better,
        reliable_change = diff_crit,
        plot = plot
    )

    output
}
## This is my ultimate testing section for maximum comfort in debugging new
## features of `rciplot()`:
# load('data/sample_data.rda')
# result <- rciplot(
#     data = sample_data,
#     pre = 'pre_data',
#     post = 'post_data',
#     reliability = 0.88,
#     recovery_cutoff = 104,
#     show_classification_percentages = T,
#     higher_is_better = T,
#     opacity = 1,
#     draw_meanmid_line = T,
#     draw_2sd_functional_line = T,
#     draw_2sd_dysfunctional_line = T,
#     mean_functional = 115,
#     mean_dysfunctional = 86,
#     sd_functional = 7.5,
#     sd_dysfunctional = 7.5
# )
# result$plot

utils::globalVariables(c(
    'category',
    'col_pre',
    'col_post',
    'intercept',
    'name',
    'slope'
))
