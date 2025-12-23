library(forestplot)
library(grid)
library(svglite)
library(rsvg)


#### Prepare labels ####
pa_labels_gf = list()
pa_labels_mt = list()

for (i in 1:nrow(plot_gf)) {
  name = plot_gf$Polyamines[i]
  if (name == 'N1-acetylspermidine') {
    pa_labels_gf[[i]] = expression(paste('N'^'1', '-acetylspermidine'))
    pa_labels_mt[[i]] = expression('')
  } else if (name == 'N8-acetylspermidine') {
    pa_labels_gf[[i]] = expression(paste('N'^'8', '-acetylspermidine'))
    pa_labels_mt[[i]] = expression('')
  } else {
    pa_labels_gf[[i]] = name
    pa_labels_mt[[i]] = expression('')
  }
}


#### Create table####
tabletext_gf = list(
  c(list(expression(bold('Polyamines'))), pa_labels_gf),
  c('HR [95% CI]', plot_gf$HR_CI),
  c('P value', plot_gf$P_fmt)
)

tabletext_mt = list(
  c(list(expression('')), pa_labels_mt),
  c('HR [95% CI]', plot_mt$HR_CI),
  c('P value', plot_mt$P_fmt)
)


#### Generate combined forest plot ####
svglite('combined_forest_plot.svg', width = 9, height = 6)

grid.newpage()
pushViewport(viewport(width = 0.92, height = 0.80, x = 0.5, y = 0.5))

# Add legend
x_pos = 0.9
y_pos = 0.85

grid.rect(x = x_pos, y = y_pos, width = unit(35, 'mm'), height = unit(12, 'mm'),
         gp = gpar(fill = 'white', col = 'black', lwd = 0.5))

grid.rect(x = x_pos - 0.06, y = y_pos + 0.015,
         width = unit(1.5, 'mm'), height = unit(1.5, 'mm'),
         gp = gpar(fill = '#009E73', col = '#009E73'))
grid.rect(x = x_pos - 0.06, y = y_pos - 0.015,
         width = unit(1.5, 'mm'), height = unit(1.5, 'mm'),
         gp = gpar(fill = '#d55e00', col = '#d55e00'))

grid.text('Graft Failure', x = x_pos - 0.03, y = y_pos + 0.015,
         just = 'left', gp = gpar(fontsize = 10, font = 2))
grid.text('Mortality', x = x_pos - 0.03, y = y_pos - 0.015,
         just = 'left', gp = gpar(fontsize = 10, font = 2))


pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(0.65, 0.35), 'npc')),
                     y = 0.43, height = 0.86, just = 'center'))

# Graft failure forest plot
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
forestplot(
  tabletext_gf,
  graph.pos = 2,
  mean = c(NA, plot_gf$HR),
  lower = c(NA, plot_gf$CI_lower),
  upper = c(NA, plot_gf$CI_upper),
  xlab = 'Hazard Ratio [95% CI]',
  txt_gp = fpTxtGp(label = gpar(cex = 0.8, fontfamily = 'sans'),
                   ticks = gpar(cex = 0.7),
                   xlab = gpar(cex = 0.8, font = 2)),
  col = fpColors(box = '#009E73', lines = '#009E73', zero = 'grey60'),
  zero = 1,
  boxsize = 0.2,
  ci.vertices = TRUE,
  ci.vertices.height = 0.1,
  lwd.ci = 1,
  grid = TRUE,
  xticks = seq(0.4, 1.8, by = 0.2),
  clip = c(0.4, 1.8),
  new_page = FALSE,
  is.summary = c(TRUE, rep(FALSE, nrow(plot_gf))),
  hrzl_lines = list('2' = gpar(lty = 2)),
  fn.ci_norm = fpDrawNormalCI,
  vertices = TRUE,
  graphwidth = unit(5, 'cm'),
  colgap = unit(0.5, 'mm'),
  cex = 0.9,
  mar = unit(c(0, 0, 0, 0), 'mm'),
  lineheight = unit(8, 'mm')
)
popViewport()

# Mortality forest plot
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
forestplot(
  tabletext_mt,
  graph.pos = 2,
  mean = c(NA, plot_mt$HR),
  lower = c(NA, plot_mt$CI_lower),
  upper = c(NA, plot_mt$CI_upper),
  xlab = 'Hazard Ratio [95% CI]',
  txt_gp = fpTxtGp(label = gpar(cex = 0.8, fontfamily = 'sans'),
                   ticks = gpar(cex = 0.7),
                   xlab = gpar(cex = 0.8, font = 2)),
  col = fpColors(box = '#d55e00', lines = '#d55e00', zero = 'grey60'),
  zero = 1,
  boxsize = 0.2,
  ci.vertices = TRUE,
  ci.vertices.height = 0.1,
  lwd.ci = 1,
  grid = TRUE,
  xticks = seq(0.4, 1.4, by = 0.2),
  clip = c(0.4, 1.4),
  new_page = FALSE,
  is.summary = c(TRUE, rep(FALSE, nrow(plot_mt))),
  hrzl_lines = list('2' = gpar(lty = 2)),
  fn.ci_norm = fpDrawNormalCI,
  vertices = TRUE,
  graphwidth = unit(4, 'cm'),
  colgap = unit(0.5, 'mm'),
  cex = 0.9,
  mar = unit(c(0, 0, 0, 8), 'mm'),
  lineheight = unit(8, 'mm')
)
popViewport(2)

dev.off()

rsvg_pdf('combined_forest_plot.svg', 'combined_forest_plot.pdf')
