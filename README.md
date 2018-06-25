


# A functional landscape of chronic kidney disease entities from public transcriptomic data


To develop efficient therapies and identify novel early biomarkers for chronic kidney disease (CKD) an understanding of the molecular mechanisms orchestrating it is essential. We here set out to understand how differences in CKD origin are reflected in gene regulatory mechanisms. To this end, we collected and integrated publicly available human-kidney glomerular microarray gene expression data for nine kidney disease entities that account for a majority of CKD worldwide [Focal segmental glomerulosclerosis (FSGS), Minimal Change Disease (MCD), FSGS-MCD, IgA nephropathy (IgAN), Lupus nephritis (LN), Membranous glomerulonephropathy (MGN), Diabetic nephropathy (DN), Hypertensive nephropathy (HN) and Rapidly progressive glomerulonephritis (RPGN)]. We included data from five distinct studies and compared glomerular gene expression profiles to that of non-tumor part of kidney cancer nephrectomy tissues. A major challenge was the integration of the data from different sources, platforms and conditions, that we mitigated with a bespoke stringent procedure. This allowed us to perform a global transcriptome-based delineation of different kidney disease entities, obtaining a landscape of their similarities and differences based on the genes that acquire a consistent differential expression between each kidney disease entity and tumor nephrectomy. Furthermore, we derived functional insights by inferring signaling pathway and transcription factor activity from the collected gene expression data, and identified potential drug candidates based on expression signature matching. These results provide a foundation to comprehend the specific molecular mechanisms underlying different kidney disease entities, that can pave the way to identify biomarkers and potential therapeutic targets.

The corresponding article for this project is available on [bioRxiv (pdf)](https://www.biorxiv.org/content/biorxiv/early/2018/02/14/265447.full.pdf).

```
@article {Tajti265447,
	author = {Tajti, Ferenc and Antoranz, Asier and Ibrahim, Mahmoud M. and Kim, Hyojin and Ceccarelli, Francesco and Kuppe, Christoph and Alexopoulos, Leonidas G. and Kramann, Rafael and Saez-Rodriguez, Julio},
	title = {A functional landscape of chronic kidney disease entities from public transcriptomic data},
	year = {2018},
	doi = {10.1101/265447},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2018/02/14/265447},
	eprint = {https://www.biorxiv.org/content/early/2018/02/14/265447.full.pdf},
	journal = {bioRxiv}
}
```

<p align="center">
    Flow of analysis 
    <img src="https://github.com/saezlab/CKD_Landscape/blob/master/Plot/Analysis_workflow.png" width="900" height="300">
</p>



## Batch effect mitigation and platform integration 

* Setting 1: **bachMitigationGPL96.R**
* Setting 2: **bachMitigationGPL570.R**
* Setting 3: **platformIntegration.R**

## Differential expression and meta analysis 

* Setting 4: **differentialAnalysis.R**

## Pathway analysis with PROGENy

* Setting 5: **progenyPathwayAnalysis.R**

## Transcription Factor activity analysis with DoRothEA

* Setting 6: **transcriptionFactorsDorothea.R**

## Diffusion Map 

* Setting 7: **DiffusionMap.R**

<p align="center">
    Similarity Landscape of CKD entities 
    <img src="https://github.com/saezlab/CKD_Landscape/blob/master/Plot/diffusion_gex_glom.png" width="700" height="700">
</p>


For colour coding refer to [bioRxiv(pdf)](https://www.biorxiv.org/content/biorxiv/early/2018/02/14/265447.full.pdf) Figure 2B

## Drug repositioning


Drug repositioning has two main part, first is to get characteristics direction and second, to find drug candidates for each CKD entity. 

* First, in order to get the characteristic direction (the signature of cosine distances) for each CKD entity 
we used the normalised gene expression values and then we measured the signature of cosine distances 
with the open source code from [Implementations of the Characteristic Direction Method](http://maayanlab.net/CD/)

* Second, to identify drug candidates from drug repositioning,
we utilized the signature of cosine distances from the first step as an input data and 
applied them to the search engine [L1000 Characteristic Direction Signature Search Engine](http://amp.pharm.mssm.edu/L1000CDS2/#/index)

The above steps were repeated for each CKD entity, separately. 


## License

Distributed under the GNU GPLv3 License. See accompanying file LICENSE.txt or copy at http://www.gnu.org/licenses/gpl-3.0.html.








