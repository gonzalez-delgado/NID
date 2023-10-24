# STINA: Statistical Testing of the Interdependence between the effects of Neighboring Amino Acids

STINA proposes a statistical methodology to test whether the influences of left and right neighbor identities on the backbone dihedral angles of a given amino acid residue are independent. In other words, STINA assesses if left and right neighbors can be separately taken into account when estimating the (phi, psi) distribution for the central residue or, on the contrary, if the identity of both neighbors has to be simultaneously considered.

The presented methodology is based on a statistical test assessing the null hypothesis "the influences of the two neighbors are independent". Small (close to zero) p-values are obtained if the given samples have very little chance of being drawn if the null hypothesis is true. This suggests that both influences are interdependent, and therefore that the null hypothesis can be rejected with statistical guarantees.

The function ``independence_test`` -included in the [independence_test.R](independence_test.R) file together with its implementation guidelines- performs the independence test given a sample of (phi, psi) values for a given central amino acid with different combinations on left and right neighbors. The file [ALA_angles.RData](ALA_angles.RData) is an example of such data, containing a sample of dihedral angles of an alanine amino acid. Some examples of analyses performed on [ALA_angles.RData](ALA_angles.RData) data are included in [examples.R](examples.R). The first corresponds to the implementation of the statistical test with three different discretization methods. The second evaluates the effect of neighbor polarity on the magnitude of interdependence, for different neighbor sets.

Details on the method and results are provided in:

J. González-Delgado, Pau Bernadó, Pierre Neuvial, Juan Cortés. Statistical proofs of the interdependence between nearest neighbor effects on polypeptide backbone conformations. <i>Journal of Structural Biology</i>, 214(4): 107907, 2022. [https://doi.org/10.1016/j.jsb.2022.107907](url).

The database used in the work is available at [https://moma.laas.fr/static/data/tripeptide_angles_data.tar](https://moma.laas.fr/static/data/tripeptide_angles_data.tar).
