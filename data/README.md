# Dataset and experimental units

`catalyst_activity.csv` contains the 24 observations supplied by Alessio Cosenza for this learning project. The final HTML line breaks in the supplied table were removed; numerical values were preserved. No observations were removed, averaged away, imputed, or recalculated in the raw file.

There are 12 catalyst samples: 2 nominal Pt loadings × 3 carbon supports × 2 Pt localities. Each sample has two **separately prepared electrode tests of the same catalyst batch**, as confirmed by the data owner. These are not independent catalyst syntheses.

| Column | Meaning | Treatment in this project |
|---|---|---|
| `Sample` | Catalyst identifier | Source of nominal loading and carbon labels |
| `Loading` | Measured Pt loading, interpreted as wt% Pt in the catalyst | Preserved as a continuous measurement; not the primary factorial predictor |
| `Locality` | IN or OUT Pt placement category | Categorical factor |
| `Replicate` | Electrode test 1 or 2 within a catalyst sample | Row identifier; not automatically a block or paired-test indicator |
| `ECSA_CO` | CO-derived electrochemically active surface area | Convention: m²/gPt; retained for audit |
| `ECSA_TEM` | TEM-derived surface area | Convention: m²/gPt; not used as a predictor |
| `MA` | Mass activity | Convention: A/gPt; retained for audit |
| `SA` | Supplied specific activity | Response; convention: µA/cm²Pt |

The supplied table did not explicitly state units. The conventions above are inferred from the electrochemical quantities and the approximate relationship `SA = 100 × MA / ECSA_CO`. Verify these conventions against laboratory records before manuscript use. The factor called loading here is the Pt mass fraction in the catalyst, not electrode areal Pt loading.

Derived columns are `Loading_level` (15 or 40, extracted from the sample identifier), `Carbon` (Disordered, Graphitized, FCX-B), and audit columns. Measured loadings span 12.2–18.8 wt% in the nominal 15 group and 38.4–45.01 wt% in the nominal 40 group. They are constant across the two tests of each catalyst.

Some supplied SA values differ slightly from `100 × MA / ECSA_CO`; the largest absolute discrepancy is about 0.434 µA/cm²Pt. We preserve the supplied response and export the differences. Rounding or upstream processing may explain discrepancies, but their cause has not been verified. The ratio audit does not validate the underlying measurements.

No synthesis batch IDs, electrode preparation dates, run order, instrument sessions, or randomization records were supplied. Replicate 1 across different samples is not assumed to mean the same day or experimental block. The model assumes independent electrode errors conditional on each tested sample and a common error variance. Batch-to-batch variation and synthesis-level causal effects cannot be estimated from this table.
