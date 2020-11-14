# geostats

Lecture notes, example datasets and code for the introductory
statistics module for geoscientists at University College London
(UCL). The notes can be found
[here](https://github.com/pvermees/geostats/blob/main/latex/geostats.pdf)

## installation

The **geostats** package requires **R**, which can be downloaded from
[CRAN](https://www.r-project.org/). In **R**, the **geostats** package
can be installed from GitHub using the **remotes** package:

```
install.packages('remotes')
remotes::install_github('pvermees/geostats/package')
```

## example

The following code snippet creates a fifth order Koch snowflake:

```
library('geostats')
koch(n=5)
```

A complete list of documented **geostats** functions can be viewed by
entering `help(package='geostats')` at the **R** command prompt.

## Author

[Pieter Vermeesch](https://www.ucl.ac.uk/~ucfbpve/)

## License

This project is licensed under the GPL-3 License
