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



cfnBinoms = (n, Q) -> (

    binoms := phyloToricFP(leafTree(5, {{0,1}, {0,1,2}}), CFNmodel);
    liftBinoms := map(Q, ring binoms, apply(gens(ring binoms), l -> q_(prepend(0, last baseName l))));

    liftBinoms(binoms) 
    )



sunletParam = n -> (

  indR := indexSet n;

  S := QQ[{t} | apply(2*n, i -> a_i)];

  --network param
  images := for ivec in indR list (

    product(apply(n, j -> if ivec_j == 1 then a_j else 1))*(product(apply(n-1, j -> if odd sum(ivec_(toList(0..j))) then a_(n+j) else 1)) + product(apply(toList(1..n-1), j -> if odd sum(ivec_(toList(1..j))) then a_(n+j) else 1)))
    );

  return t*images;
  )

-- Creates the parameterization of the sunlet network variety in the Fourier coordinates
-- Note that we give each q_g variable degree 2n so that the elimination ideal is properly homogenized
-- This greatly speeds up the Grobner basis computation and also allows us to compute a degree-bounded Grobner basis with degLimit option
-- Note that degLimit should just be the degree that one wants the resulting polynomials in the q_g to be if each q_g has degree 1
sunletElimIdeal = {degLimit => null, Qring => null} >> opts -> n -> (

  indR := toList delete(null, apply((n:0)..(n:1), i ->  if sum(toList i) % 2 == 0 then i));

  S := QQ[{t} | apply(2*n, i -> a_i) | apply(indR, i -> q_i), Degrees =>  toList(2*n+1 : 1)|toList(2^(n-1) : 2*n - 1), MonomialOrder => {2*n+1, 2^(n-1)}];

  --network param
  images := for ivec in indR list (

    product(apply(n, j -> if ivec_j == 1 then a_j else t))*(product(apply(n-1, j -> if odd sum(ivec_(toList(0..j))) then a_(n+j) else t)) + product(apply(toList(1..n-1), j -> if odd sum(ivec_(toList(1..j))) then a_(n+j) else t)))
    );

  elimIdeal := ideal(for i from 0 to #images-1 list q_(indR#i) - images#i);
  G := selectInSubring(1, if opts.degLimit === null then gens gb(elimIdeal) else gens gb(elimIdeal, DegreeLimit => 2*n*(opts.degLimit)));
  R := if opts.Qring === null then qRing n else opts.Qring;
  return(sub(ideal G, R));
  )


n = 6
Q = qRing(n, CFNmodel);
K = sunletElimIdeal(6, degLimit => 2, Qring => Q);

X = xRing n;
psi = pfaffianMap(n, Q, X)
kerPsi = ker(psi)
I = pfaffianDetIdeal(n, X)
qCoordsInOrder = for g in gens(X) list(

    ind = last baseName(g);

    basicQCoord(n, ind#0, ind#1, Q)
    );

invPsi = map(Q, X, qCoordsInOrder)

J1 = homogenize(kerPsi, Q_0)
J2 = invPsi(I)
J = J1 + J2


for g in K_* list if g % J == 0 then continue else g




