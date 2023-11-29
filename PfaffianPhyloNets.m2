-- This file contains functions for creating the polynomial maps which parameterize the CFN sunlet network variety in different coordinate systems
-- 
restart

needsPackage "gfanInterface"
needsPackage "PhylogeneticTrees"


-- n, the number of leaves
-- i, j the positions of the two ones
-- Q, the ring of Fourier coordinates
-- Returns the Fourier coordinatr q_(e_i + e_j) where the indexing corresponds to leaves labelled 1 to n
basicQCoord = (n, i, j, Q) -> sub(q_(apply(0..(n-1), l -> if member(l, {i-1, j-1}) then 1 else 0)), Q)


-- n, the number of leaves
-- Returns the ring X_n in whose variables correspond to the entries of a generic antisymmetric matrix
-- We construct this ring with the term order described in Section 3, directly above Lemma 3.7
xRing = n -> (

    firstRowInds := apply(toList(2..n), j -> (1, j));

    otherInds := flatten reverse for i from 2 to n list for j from i+1 to n list (i, j);


    QQ[apply(firstRowInds, s -> x_s) | apply(otherInds, s -> x_s), MonomialOrder => Lex]
    );


-- n, the number of leaves
-- X, the ring X_n whose variables are of the form x_(i, j) for i > j
-- Returns a generic antisymmetric matrix Omega whose entries are given by the variables in X
omegaMatrix = (n, X) -> (

    Omega := for i from 1 to n list(
                for j from 1 to n list(

                    if j <= i then 0 else sub(x_(i, j), X)
                    )

                );

    Omega = matrix(Omega);

    Omega - transpose(Omega)
    )

-- n, the number of leaves
-- Q, the ring Q_n of Fourier coordinates
-- X, the ring X_n whose variables are of the form x_(i, j) for i > j
-- Returns the map psi_n which sends the coordinate q_g to 1/2^k Pf(Omega_g^n)
pfaffianMap = (n, Q, X) -> (

    Omega := omegaMatrix(n, X);

    psi := for qg in gens(Q) list(

            g := last baseName qg;

            S := delete(null, apply(toList(0..n-1), i -> if g_i == 1 then i));

            pf := if #S > 1 then (radical ideal det(Omega_S^S))_0 else 1
            );

    return map(X, Q, psi)
    )


-- n, the number of leaves
-- X, the ring X_n whose variables are of the form x_(i, j) for i > j
-- Returns the determinantal ideal generated by the set G_n of 2x2 and 3x3 minors described in Theorem 3.5
pfaffianDetIdeal = (n, X) -> (

    Omega := omegaMatrix(n, X);

    row := subsets(toList(1..(n - 1)), 2);

    minors2 := for r in row list(

        col := subsets((r_1 + 1)..(n - 1), 2);

        if #col == 0 then break;

        for c in col list(

            minors(2, Omega_c^r)
            )
        );

    minors2 = sum flatten minors2;

    if n < 6 then return minors2;

    minors3 := for r in row list(

        col := subsets((r_0 + 1)..(r_1 - 1), 3);

        if #col == 0 then continue;

        for c in col list(

            minors(3, Omega_c^({0}|r))
            )
        );

    minors3 = sum flatten minors3;

    return (minors2 + minors3);
    )


