# GrowingSOM
Implementation of Growing Self Organizing Maps in R

Mainly based on:
__Damminda Alahakoon, Saman K. Halgamuge (2000)__: Dynamic Self-Organizing Maps with Controlled Growth for Knowledge Discovery. IEEE TRANSACTIONS ON NEURAL NETWORKS, VOL. 11.

## Functionality

- train.gsom()
- train_xy.gsom()
- map.gsom()
- predict.gsom()
- summary.gsom()
- print.gsom()
- plot.gsom()

## Installation

Make sure you are in the package directory then run:

```bash
	R CMD build .
```

Please note that dependencies may need to be installed in advance (`install.packages("plotrix")`)

Windows only: A binary Package needs to be obtained by: 

```bash
	mkdir temp_location
	R CMD INSTALL --build -l temp_location GrowingSOM_VERSION.tar.gz
```

installation Linux (in R):
```R
	install.packages(repos=NULL,"PHATH_TO_HERE/GrowingSOM_VERSION.tar.gz")
```

installation Windows (in R):
```R
	install.packages(repos=NULL, "PATH_TO_HERE\GrowingSOM_VERSION.zip", type="binary")
```

## Legal

GrowingSOM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

GrowingSOM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GrowingSOM.  If not, see <http://www.gnu.org/licenses/>.
