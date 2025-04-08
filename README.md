# Modlling and Simulation project

_Universite Paris-Saclay, 2025_
This repository contains the code and documentation for the Modelling and Simulation project at Université Paris-Saclay. The project focuses on the development of mathematical models in biology, in here replicating a paper that have modelled tumor growth and immune interaction using agent-based modelling.

> "An agent-based model of avascular tumor growth:
> Immune response tendency to prevent cancer development" [^1]

| Tumor growth | Tumor-immune interaction |
|:--:|:--:|
| ![g](./assets/growth.gif) | ![i](./assets/immune.gif) |

[^1]: Pourhasanzade, F., Sabzpoushan, S. H., Alizadeh, A. M., & Esmati, E. (2017). An agent-based model of avascular tumor growth: Immune response tendency to prevent cancer development. Simulation, 93(8), 641-657.

You'll find the report [in docs folder](./docs/JoelleASSY_RayaneADAM_MS_report.pdf).

Models are implemented in NetLogo and found in the `models` folder. The models are:
- `tumor_growth.nlogo`: A model of tumor growth
- `tumor_immune.nlogo`: A model of tumor-immune interaction
- `tumor_growth_contribution.nlogo`: small contribution added
- `tumor_immune_contribution.nlogo`: small contribution added

