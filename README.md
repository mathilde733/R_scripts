## Epidemiological & Biological Modeling in R

This repository contains research and coursework scripts focused on mathematical modeling in biology, split into two distinct methodological approaches.

### 1. Computational Epidemiology

Script: CompEpi_HumanMosquitoDynamics.R \
Collaboration with Janosch Imhof, Florence Kurz, Lara Land and Valentin Müller. 

This module explores modern computational approaches to disease spread, specifically focusing on [vector-borne diseases](https://www.who.int/news-room/fact-sheets/detail/vector-borne- diseases) (e.g., Malaria, Dengue).

#### Model Architecture

The baseline mathematical model consists of an SIR model for human population and SI model for vector population. The model equations are adapted from a previous paper written Tang et al., 2023. \

Inspiration Schema: \
<img src="https://andreashandel.github.io/IDEMAbook/images/vectorborne-schematic.png" alt="Schema" width="500">

This dual model contains the following underlying assumptions: mosquitoes (vector, 𝑣) can be either susceptible (𝑆) or infected (𝐼) and have a short lifespan, while humans (host, h) can be either susceptible (𝑆), infected (𝐼) or recovered (𝑅). Additionally, and importantly, there is no intra-species transmission.

#### Key Features:

- Dynamic interaction between host and vector populations.

- Simulation of intervention measures: Vaccination, insecticide-treated nets (ITNs), and vector birth/death rate controls. 

#### Tools used: deSolve for numerical integration and ggplot2 for high-fidelity time-series visualization.

#### References
1. **Pagendam, D. E. et al.** (2020). [Modelling the Wolbachia incompatible insect technique: strategies for effective mosquito population elimination](https://doi.org/10.1186/s12915-020-00891-5). *BMC Biol.* 18, 161.
2. **Bhatt, S. et al.** (2013). [The global distribution and burden of dengue](https://doi.org/10.1038/nature12060). *Nature* 496, 504–507.
3. **Chen, L. H. et al.** (2023). [Epidemiology and burden of dengue fever in the United States: a systematic review](https://doi.org/10.1093/jtm/taad127). *J. Travel Med.* 30, taad127.
4. **Ioos, S. et al.** (2014). [Current Zika virus epidemiology and recent epidemics](https://doi.org/10.1016/j.medmal.2014.04.008). *Médecine Mal. Infect.* 44, 302–307.
5. **Hills, S. L., Fischer, M. & Petersen, L. R.** (2017). [Epidemiology of Zika Virus Infection](https://doi.org/10.1093/infdis/jix434). *J. Infect. Dis.* 216, S868–S874.
6. **Ribeiro Dos Santos, G. et al.** (2025). [Global burden of chikungunya virus infections and the potential benefit of vaccination campaigns](https://doi.org/10.1038/s41591-024-03221-8). *Nat. Med.* 31, 2342–2349.
7. **Tang, T.-Q. et al.** (2023). [Analysis of the dynamics of a vector-borne infection with the effect of imperfect vaccination from a fractional perspective](https://doi.org/10.1038/s41598-023-41310-7). *Sci. Rep.* 13, 14398.
8. **Handel, A.** [Infectious Disease Epidemiology: A Model-based Approach (IDEMA)](https://andreashandel.github.io/IDEMAbook/vector-borne-transmission-1.html).
9. **Rivero, A. et al.** (2010). [Insecticide Control of Vector-Borne Diseases: When Is Insecticide Resistance a Problem?](https://doi.org/10.1371/journal.ppat.1001000) *PLoS Pathog.* 6, e1001000.

### 2. Classical Models in Biology

Script: report_classical_modeling.R
Collaboration with Florence Kurz.

This module covers fundamental ecological interactions and the mathematical foundations of population dynamics.

<img src="https://andreashandel.github.io/IDEMAbook/images/vectorborne-schematic.png" alt="Schema" width="500">

#### Topics Covered:

- Lotka-Volterra Competition: Analysis of stable coexistence, competitive exclusion, and priority effects (Alternative Stable States).

- Lotka-Volterra Mutualism: Modeling positive species interactions and determining the conditions for stable vs. unstable (runaway) equilibrium points.

Mathematical Methods: 
- Phase space analysis and vector field visualization
- Nullclines (Zero-growth isoclines) calculation to determine population stability.

#### Constraints: Developed strictly using the deSolve library for core numerical evaluations.

#### References
1.  **Volterra, V.** (1926). [Fluctuations in the Abundance of a Species considered Mathematically](https://doi.org/10.1038/118558a0). *Nature*, 118, 558–560.
2.  **Hardin, G.** (1960). [The Competitive Exclusion Principle](https://doi.org/10.1126/science.131.3409.1292). *Science*, 131(3409), 1292–1297.
3.  **Chesson, P.** (2000). [Mechanisms of Maintenance of Species Diversity](https://doi.org/10.1146/annurev.ecolsys.31.1.343). *Annual Review of Ecology and Systematics*, 31, 343–366.
4.  **Fukami, T.** (2015). [Historical Contingency in Community Assembly: Integrating Niches, Priority Effects, and Species Pools](https://doi.org/10.1146/annurev-ecolsys-110411-160340). *Annual Review of Ecology, Evolution, and Systematics*, 46, 1–23.
5.  **Boucher, D. H., James, S., & Keeler, K. H.** (1982). [The Ecology of Mutualism](https://doi.org/10.1146/annurev.es.13.110182.001531). *Annual Review of Ecology and Systematics*, 13, 315–347.
6.  **Stevens, M. H. H.** (2009). [A Primer of Ecology with R](https://doi.org/10.1007/978-0-387-89882-7). *Springer*.
