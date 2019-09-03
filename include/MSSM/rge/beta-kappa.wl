{
 {
  "kappa", 
  (
  - 6 (g1^2 + 5 g2^2 - 5 trace[Yu, Adj[Yu]]) kappa 
  - 6 (g1^2 + 5 g2^2 - 5 trace[Yu, Adj[Yu]]) Tp[kappa] 
  + 5 (
        MatMul[kappa, conj[Ye], Tp[Ye]] 
      + MatMul[Ye, Adj[Ye], kappa] 
      + MatMul[Ye, Adj[Ye], Tp[kappa]] 
      + MatMul[Tp[kappa], conj[Ye], Tp[Ye]]
      )
  ) / 10, 
  (
    (
      40 (g1^2 + 20 g3^2) trace[Yu, Adj[Yu]] 
    + 3 (
          69 g1^4 
        + 30 g1^2 g2^2 
        + 125 g2^4 
        - 50 trace[Yd, Adj[Yd], Yu, Adj[Yu]] 
        - 150 trace[Yu, Adj[Yu], Yu, Adj[Yu]]
        )
    ) kappa 
  + (
      40 (g1^2 + 20 g3^2) trace[Yu, Adj[Yu]] 
    + 3 (
          69 g1^4 
        + 30 g1^2 g2^2 
        + 125 g2^4 
        - 50 trace[Yd, Adj[Yd], Yu, Adj[Yu]] 
        - 150 trace[Yu, Adj[Yu], Yu, Adj[Yu]]
        )
    ) Tp[kappa] 
  + 5 (
        (6 g1^2 - 15 trace[Yd, Adj[Yd]] - 5 trace[Ye, Adj[Ye]]) MatMul[kappa, conj[Ye], Tp[Ye]] 
      + (6 g1^2 - 15 trace[Yd, Adj[Yd]] - 5 trace[Ye, Adj[Ye]]) MatMul[Ye, Adj[Ye], kappa] 
      + 6 g1^2 MatMul[Ye, Adj[Ye], Tp[kappa]] 
      - 15 trace[Yd, Adj[Yd]] MatMul[Ye, Adj[Ye], Tp[kappa]] 
      - 5 trace[Ye, Adj[Ye]] MatMul[Ye, Adj[Ye], Tp[kappa]] 
      + 6 g1^2 MatMul[Tp[kappa], conj[Ye], Tp[Ye]] 
      - 15 trace[Yd, Adj[Yd]] MatMul[Tp[kappa], conj[Ye], Tp[Ye]] 
      - 5 trace[Ye, Adj[Ye]] MatMul[Tp[kappa], conj[Ye], Tp[Ye]] 
      - 10 MatMul[kappa, conj[Ye], Tp[Ye], conj[Ye], Tp[Ye]] 
      - 10 MatMul[Ye, Adj[Ye], Ye, Adj[Ye], kappa] 
      - 10 MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Tp[kappa]] 
      - 10 MatMul[Tp[kappa], conj[Ye], Tp[Ye], conj[Ye], Tp[Ye]]
      )
  ) / 50
 }
}
