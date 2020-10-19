# network-significance

Script to estimate P-values of network metrics and compare pairs of networks using Monte Carlo procedures in R.

[Ecological Synthesis Lab](https://marcomellolab.wordpress.com) (SintECO).

Authors: Renata Muylaert, Pavel Dodonov & Marco Mello.

E-mail: renatamuy@gmail.com.

First published on April 25th, 2017 (English version).

Run in R version 4.0.2 (2020-06-22) -- "Taking Off Again".

<a href="https://doi.org/10.5281/zenodo.1487607"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1487607.svg" alt="DOI"></a>

Disclaimer: You may use this software freely for any purposes at your own risk. We assume no responsibility or liability for the use of this software, convey no license or title under any patent, copyright, or mask work right to the product. We reserve the right to make changes in the software without notification. We also make no representation or warranty that such application will be suitable for the specified use without further testing or modification. If this software helps you produce any academic work (paper, book, chapter, dissertation, report, etc.), please acknowledge the authors and cite this repo's URL.


## Functionality and origin

This script provided a series of codes to estimate the significance (P-value) of network-level metrics for a single network and for comparing two networks. We first wrote it to use in a series of papers published by our lab since 2016, and have been updating it since then to make the code shorter and more efficient.


## List of files

1. net1.txt -> example data from Bezerra et al. (2009, J.Anim.Ecol.). Interactions of pollination between oil-collecting bees (Centridini) and oil flowers (Malpighiaceae).

2. net2.txt -> example data from Queiroz et al. (2020, Biotropica). Interactions of pollination between bats (Chiroptera), hawkmoths (Sphingidae), and plants.

3. network-significance.R -> script to estimate P-values of network metrics in different situations. This script includes code to plot the distribution of randomized values and the observed value. Always visualize the data before calculating the statistics!


## Instructions

Follow the instructions provided in the script "network-significance.R".


## Feedback

If you have any questions, suggestions, or corrections, please feel free to open an [issue](https://github.com/marmello77/network-significance/issues) or make a [pull request](https://github.com/marmello77/network-significance/pulls).


## Data sources

1. net1 -> Bezerra, E. L. S., Machado, I. C., & Mello, M. A. R. (2009). [Pollination networks of oil-flowers: a tiny world within the smallest of all worlds](https://doi.org/10.1111/j.1365-2656.2009.01567.x). Journal of Animal Ecology, 78(5), 1096–1101. 

2. net2 -> Queiroz, J. A., Diniz, U. M., Vázquez, D. P., Quirino, Z. M., Santos, F. A. R., Mello, M. A. R., & Machado, I. C. (2020). Bats and hawkmoths form mixed modules with flowering plants in a nocturnal interaction network. Biotropica, *accepted*. See also this [repo](https://github.com/marmello77/queiroz-et-al-2020).


## Acknowledgments

We thank our labmates and our sponsors, especially the Alexander von Humboldt-Stiftung, CNPq, CAPES, and FAPESP, who gave us grants, fellowships, and scholarships. Last, but not least, we thank the [Stack Overflow Community](https://stackoverflow.com), where we solve most of our coding dilemmas. 


## Suggested readings

If you want to understand the concepts used in our script, read the following works:

* Beckett, S. J. (2016). [Improved community detection in weighted bipartite networks](https://doi.org/10.1098/rsos.140536). Royal Society Open Science, 3(1), 140536.

* Blüthgen, N., Menzel, F., Hovestadt, T., Fiala, B., & Bluthgen, N. (2007). [Specialization, constraints, and conflicting interests in mutualistic networks](https://doi.org/10.1016/j.cub.2006.12.039). Current Biology, 17(4), 341–346.

* Dormann, C. F., Gruber, B., & Fründ, J. (2008). [Introducing the bipartite package: analyzing ecological networks](https://www.uni-goettingen.de/de/document/download/96729eb9d30a6f2dc4403df15854305c.pdf/Rnews2008,8_8-11_open.pdf). R News, 8(2), 8–11.

* Felix, G. M., Pinheiro, R. B. P., Poulin, R., Krasnov, B. R., & Mello, M. A. R. (2017). [The compound topology of a continent-wide interaction network explained by an integrative hypothesis of specialization](https://doi.org/10.1101/236687). BioRxiv, 236687.

* Mello, M. A. R., Felix, G. M., Pinheiro, R. B. P., Muylaert, R. L., Geiselman, C., Santana, S. E., … Stevens, R. D. (2019). [Insights into the assembly rules of a continent-wide multilayer network](https://doi.org/10.1038/s41559-019-1002-3). Nature Ecology & Evolution, 3(11), 1525–1532.

* Pinheiro, R. B. P., Felix, G. M. F., Dormann, C. F., & Mello, M. A. R. (2019). [A new model explaining the origin of different topologies in interaction networks](https://doi.org/10.1002/ecy.2796). Ecology, 100(9), e02796.

* Ulrich, W., Almeida-Neto, M., & Gotelli, N. J. (2009). [A consumer’s guide to nestedness analysis](https://doi.org/10.1111/j.1600-0706.2008.17053.x). Oikos, 118(1), 3–17.
