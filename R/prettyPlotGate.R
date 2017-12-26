prettyPlotGate <- function(gs, CellTags, population = "nonMaxF") {

if (!exists("template")) template <- list("G.nRange" = c(0.5,4.5), "R.nRange" = c(0.4,4.5))


overlay <- CellTags
nOverlay <- length(overlay)
#pal <- c(brewer.pal(9,"Set1"), brewer.pal(12,"Set3"))
#CT color alphabet
pal <- c(
"#191919","#f0a3ff","#0075dc","#4c005c","#2bce48","#808080","#ffa405","#ff0010","#5ef1f2","#ff5005","#993f00","#005c31","#ffcc99","#8f7c00","#9dcc00","#c20088","#003380","#ffa8bb","#00998f","#740aff"
)
overlay.fill <- pal[1:nOverlay]
names(overlay.fill) <- overlay
overlay.symbol <- sapply(overlay.fill, function(col)list(fill = col), simplify = FALSE)

key = list(text = list(names(overlay.symbol), cex = 4)
                   , points = list(col = overlay.fill
                                   , pch = 19
                                   , cex = 4)
             , columns = 1
             , between = 0.5
             , space = "right"
             , padding.text = 5)

p <- plotGate(gs, population, strip=T, stats=F, xbin=128
		, par.settings = list(
			  gate.text = list(background = list(fill = "transparent"))
			, panel.background = list(col = "white")
			, par.xlab.text = list(cex = 2)
			, par.ylab.text = list(cex = 2)
			, axis.text = list(cex = 2)
			)
		, arrange.main = F
		, gpar = list(nrow = 3)
		, overlay = overlay
		, overlay.symbol = overlay.symbol
		, par.strip.text = list(cex = 2)
		, xlim=template$G.nRange
		, ylim=template$R.nRange
		, key = key
		)
return (p)
}
