# Large carnivore distribution maps for Europe 2017 – 2022/23

[https://doi.org/10.5061/dryad.3xsj3txrc](https://doi.org/10.5061/dryad.3xsj3txrc)

## Description of the data and file structure

The mapping approach generally follows the methods described in (Chapron et al. 2014) and (Kaczensky et al. 2013). It updates the published Species Online Layers 2012-2016 for brown bear, Eurasian lynx, wolf, golden jackal, and wolverine  (Kaczensky et al. 2021; Ranc et al. 2022) for the period 2017-2022/23.

Large carnivore presence was mapped at a 10 x 10 km (ETRS89-LAEA Europe) grid scale. This grid is widely used for Habitat Directive reporting to the European Union (EU) and can be downloaded at: [http://www.eea.europa.eu/data-and-maps/data/eea-reference-grids-2](http://www.eea.europa.eu/data-and-maps/data/eea-reference-grids-2). The map encompasses the continental EU countries plus Switzerland and Norway, and the EU candidate / potential candidate countries in the Balkan region, in addition to Ukraine and Turkey. For the two latter countries, only parts were included; for Ukraine only the Carpathian region (for this report Ukraine was artificially cut off and the straight line in the east does not represent the national border), and the European part of Turkey (Fig. 1).

For the 2012-2016 mapping, several countries were not or not fully (not for all species) included (Hungary, Montenegro, Turkey), so that no comparisons can be made of the updated carnivore distributions with those from the last mapping for these countries.

Mapping large carnivores for this report had a two-fold goal:

·         Visualizing areas of large carnivore presence

·         Visualizing the variation in the underlying data quality

### Files and variables

The files go together with the report: Large carnivore distribution maps and population updates 2017 – 2022/23. 

Kaczensky, P., Ranc, N., Hatlauf, J., Payne, J.C. *et al.* 2024. Large carnivore distribution maps and population updates 2017 – 2022/23. Report to the European Comission under contract N° 09.0201/2023/907799/SER/ENV.D.3 “Support for Coexistence with Large Carnivores”, “B.4 Update of the distribution maps”. IUCN/SSC Large Carnivore Initiative for Europe (LCIE) and Istituto di Ecologia Applicata (IEA).

The shapefiles contain the 10x10 cells occupied by brown bear (*Ursus arctos*), Eurasian lynx (*Lynx lynx*), wolf (*Canis lupus*), golden jackal (*Canis aureus*), and wolverine (*Gulo gulo*) at the European scale. The metadata of the shape files contains the following information:

| **Metadata table** | **Information provided**                                                                           |
| :----------------- | :------------------------------------------------------------------------------------------------- |
| FID                | Unique identifier ID                                                                               |
| CELLCODE           | 10x10 km ETRS89-LAEA (Lambert Azimuthal Equal Area) Europe grid ID                                 |
| EOFORIGIN          | East coordinate in ETRS89-LAEA projection ([EPSG:3035](https://epsg.io/3035))                      |
| NOFORIGIN          | North coordinate in ETRS89-LAEA projection ([EPSG:3035](https://epsg.io/3035))                     |
| COUNTRY            | Country (in some cases large transboundary region)                                                 |
| PERSON             | Main person(s) who compiled and/or sent the map                                                    |
| SPECIES            | *Canis lupus*, *Canis aureus*, *Gulo gulo*, *Lynx lynx*, or *Ursus arctos*                         |
| POPULATION         | Species-specific population as defined by LCIE                                                     |
| PRESENCE           | Presence category: Undefined, Permanent, or Sporadic                                               |
| DATAQUAL           | Data quality categories: see report Kaczensky et al. 2024                                          |
| DATASOURCE         | Short reference of data source – for details see  Kaczensky et al. 2024                            |
| YEAR               | Time period the data layer covers                                                                  |
| YRCOMPILED         | Year the maps were compiled: 2024                                                                  |
| COMPILERS          | Kaczensky, Ranc, Hatlauf, Payne *et al.* 2024 for the Large Carnivore Initiative for Europe (LCIE) |

#### File: Brown\_bear.zip

**Description:** Shape files: Distribution of the brown bear in Europe 2017-2022/23.

#### File: Eurasian\_lynx.zip

**Description:** Shape files: Distribution of the Eurasian lynx in Europe 2017-2022/23.

#### File: Golden\_jackal.zip

**Description:** Shape files: Distribution of the Golden jackal in Europe 2017-2022/23.

#### File: Wolf.zip

**Description:** Shape files: Distribution of the wolf in Europe 2017-2022/23.

#### File: Wolverine.zip

**Description:** Shape files: Distribution of the wolverine in Europe 2017-2022/23.

## Code/software

The shape files can be visualised with any Geographic Information System (GIS) software, such as ArcGIS, QGIS, and others.

## Access information

Other publicly accessible locations of the data:

* [https://www.lcie.org/Large-carnivores](https://www.lcie.org/Large-carnivores) 
* [https://environment.ec.europa.eu/topics/nature-and-biodiversity/habitats-directive/large-carnivores_en](https://environment.ec.europa.eu/topics/nature-and-biodiversity/habitats-directive/large-carnivores_en)

