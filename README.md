README
================
Andrew Manderson

rjmonopoly
==========

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
