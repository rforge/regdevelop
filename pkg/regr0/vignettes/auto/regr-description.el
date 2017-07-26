(TeX-add-style-hook
 "regr-description"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("geometry" "a4paper" "text={14.5cm,22cm}")))
   (TeX-run-style-hooks
    "latex2e"
    "regr-description-concordance"
    "article"
    "art11"
    "graphicx"
    "Sweave"
    "inputenc"
    "geometry"
    "color"
    "booktabs"
    "amsmath"
    "texab")
   (TeX-add-symbols
    '("code" 1)
    '("wb" 1)
    '("sups" 1)
    '("mx" 1)
    '("bmath" 1)
    '("vc" 1)
    '("Vneed" 1)
    '("Hneed" 1)
    '("Tit" 1)
    "T"
    "ul")
   (LaTeX-add-labels
    "fig:taordered"
    "fig:plmatrix"
    "fig:mbox"
    "fig:mboxes2"
    "app.regr")))

