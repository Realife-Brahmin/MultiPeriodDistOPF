# Reading the frozen inertia diagnostic correctly

The replay records eigenvalues of K_s=(K+K')/2, whereas production Cholesky
uses Symmetric(K), i.e. the matrix K_u completed from K's upper triangle.
These are identical for exact symmetry, not necessarily bit-identical for
the raw late-iterate arrays. Do not silently equate their numerical inertias.

The collector adds a conservative certificate without another replay. Let
D=K-K', alpha=norm(D,F)/max(1,norm(K,F)), a=max(abs(eig(K_s))), and n=size(K,1).
The diagonal is unchanged and norm(K_u-K_s,F)=norm(D,F)/2. Also symmetric and
antisymmetric parts are Frobenius-orthogonal. For alpha<2 this implies

    norm(K_u-K_s,2) <= alpha/2 * max(1, sqrt(n)*a/sqrt(1-alpha^2/4)).

If the smallest absolute recorded eigenvalue (a/absolute_condition) exceeds
this bound plus the recorded numerical-zero threshold, perturbation cannot
change any eigenvalue sign within that margin. The reported
K_s inertia then also certifies K_u's mathematical signs, subject to ordinary
floating-point eigenvalue accuracy. This is a conservative diagnostic, not
an interval-arithmetic proof or a certificate of optimizer convergence.
If the bound fails or conditioning is nonfinite, retain an UNCERTIFIED record;
do not change data, tolerance, or declare production inertia preserved.

triangle_certificates.csv contains every bound/gap/result. The summary counts
uncertified records per mode/point. Original raw inertia, asymmetry, extrema
and condition estimates remain unchanged in the archived replay CSVs.
