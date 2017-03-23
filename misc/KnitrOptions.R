

.libPaths(c('~/R/x86_64-pc-linux-gnu-library/', .libPaths()))

# Default figure width
propOfTextwidth = 0.8

opts_chunk$set(
  echo = FALSE, 
  error = FALSE,
  cache = TRUE, 
  dev = c('postscript', 'pdf', 'png'),
  warning = FALSE,
  results = 'hide',
  message = FALSE,
  fig.width = 6.45 * propOfTextwidth,
  fig.height = 6.45 * propOfTextwidth * 0.5,
  fig.align = 'center',
  fig.pos = 't',
  #fig.showtext = TRUE, # Embeds fonts.
  out.width = paste0(propOfTextwidth,'\\textwidth'),
  dev.args=list(bg="transparent"),
  size = 'footnotesize',
  cache.extra = list(R.version)
  )
options(digits = 2)




options('scipen' = -1)

  
