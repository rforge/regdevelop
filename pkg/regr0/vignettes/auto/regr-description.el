(TeX-add-style-hook
 "regr-description"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("geometry" "a4paper" "text={14.5cm,22cm}")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "graphicx"
    "inputenc"
    "geometry"
    "color"
    "booktabs"
    "amsmath"
    "regr-desc")
   (LaTeX-add-labels
    "fig:taordered"
    "fig:plmatrix"
    "fig:mboxes"
    "fig:mboxes2")))

