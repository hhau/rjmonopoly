rjmonopoly: Automatic(ish) monotonic polynomial regression
================
Andrew Manderson

Overview
========

This package enables users to fit monotonic polynomials of varying degree using a reversible jump sampler.

Installation
------------

This package is not in CRAN right now. Just install using devtools and github

    devtools::install_github('hhau/rjmonopoly')

Overview
--------

To fit, simply call the rjmonopoly function:

    fit <- rjmonopoly::rjmonopoly(x, y, d_min = 2, d_max = 10)

The help file should be reasonably self-explanatory for the other parameters. At the moment there are only rudamentary plotting functions avaliable.

    plotDegreePost(fit)
    plotFit(fit)

I intend on adding more / some sensible generics, once I figure out what those should be.

NB
==

This sampler (and code) is currently very very delicate, and is likely to break, hang, and emit no errors. Sorry about that, I'll try to fix that in the future.

To do:
======

-   A method for finding a good starting value for the variance parameter
-   Diagnostics to identify samples that haven't performed appropriately
-   Adaptive innovation variances?
-   Multiple chains over multiple cores
-   A whole lot of these functions should be ported to `c++`, because some chains need to run for millions of iterations, and no one wants to wait for `R` to do that.
-   Flat prior option.
