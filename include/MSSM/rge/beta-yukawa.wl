{
 {
  "Ye",
    (
    - (9 g1^2) / 5 
    - 3 g2^2 
    + 3 trace[Yd, Adj[Yd]] 
    + trace[Ye, Adj[Ye]]
    ) Ye 
  + 3 MatMul[Ye, Adj[Ye], Ye] 
  + MatMul[Yv, Adj[Yv], Ye], 
    (
      (
        135 g1^4 
      + 18 g1^2 g2^2 
      + 75 g2^4 
      - 4 (g1^2 - 40 g3^2) trace[Yd, Adj[Yd]] 
      + 12 g1^2 trace[Ye, Adj[Ye]] 
      - 90 trace[Yd, Adj[Yd], Yd, Adj[Yd]] 
      - 30 trace[Yd, Adj[Yd], Yu, Adj[Yu]] 
      - 30 trace[Ye, Adj[Ye], Ye, Adj[Ye]] 
      - 10 trace[Ye, Adj[Ye], Yv, Adj[Yv]]
      ) Ye
    ) / 10 
  + (6 g2^2 - 9 trace[Yd, Adj[Yd]] - 3 trace[Ye, Adj[Ye]]) MatMul[Ye, Adj[Ye], Ye] 
  - 3 trace[Yu, Adj[Yu]] MatMul[Yv, Adj[Yv], Ye] 
  - trace[Yv, Adj[Yv]] MatMul[Yv, Adj[Yv], Ye] 
  - 4 MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye] 
  - 2 MatMul[Ye, Adj[Ye], Yv, Adj[Yv], Ye] 
  - 2 MatMul[Yv, Adj[Yv], Yv, Adj[Yv], Ye]
 }, 
 {
  "Yv", 
    (
    - (3 g1^2) / 5 
    - 3 g2^2 
    + 3 trace[Yu, Adj[Yu]] 
    + trace[Yv, Adj[Yv]]
    ) Yv 
  + MatMul[Ye, Adj[Ye], Yv] 
  + 3 MatMul[Yv, Adj[Yv], Yv], 
    (
      (207 g1^4) / 50 
    + (9 g1^2 g2^2) / 5 
    + (15 g2^4) / 2 
    + (4 (g1^2 + 20 g3^2) trace[Yu, Adj[Yu]]) / 5 
    - 3 trace[Yd, Adj[Yd], Yu, Adj[Yu]] 
    - trace[Ye, Adj[Ye], Yv, Adj[Yv]] 
    - 9 trace[Yu, Adj[Yu], Yu, Adj[Yu]] 
    - 3 trace[Yv, Adj[Yv], Yv, Adj[Yv]]
    ) Yv 
  + ((6 g1^2) / 5 - 3 trace[Yd, Adj[Yd]] - trace[Ye, Adj[Ye]]) MatMul[Ye, Adj[Ye], Yv] 
  + (6 g1^2 MatMul[Yv, Adj[Yv], Yv]) / 5 
  + 6 g2^2 MatMul[Yv, Adj[Yv], Yv] 
  - 9 trace[Yu, Adj[Yu]] MatMul[Yv, Adj[Yv], Yv] 
  - 3 trace[Yv, Adj[Yv]] MatMul[Yv, Adj[Yv], Yv] 
  - 2 MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Yv] 
  - 2 MatMul[Yv, Adj[Yv], Ye, Adj[Ye], Yv] 
  - 4 MatMul[Yv, Adj[Yv], Yv, Adj[Yv], Yv]
 }, 
 {
  "Yd", 
    (
    - (7 g1^2) / 15 
    - 3 g2^2 
    - (16 g3^2) / 3 
    + 3 trace[Yd, Adj[Yd]] 
    + trace[Ye, Adj[Ye]]
    ) Yd 
  + 3 MatMul[Yd, Adj[Yd], Yd] 
  + MatMul[Yu, Adj[Yu], Yd], 
    (
      (287 g1^4) / 90 
    + g1^2 g2^2 
    + (15 g2^4) / 2 
    + (8 g1^2 g3^2) / 9 
    + 8 g2^2 g3^2 
    - (16 g3^4) / 9 
    - (2 (g1^2 - 40 g3^2) trace[Yd, Adj[Yd]]) / 5 
    + (6 g1^2 trace[Ye, Adj[Ye]]) / 5 
    - 9 trace[Yd, Adj[Yd], Yd, Adj[Yd]] 
    - 3 trace[Yd, Adj[Yd], Yu, Adj[Yu]] 
    - 3 trace[Ye, Adj[Ye], Ye, Adj[Ye]] 
    - trace[Ye, Adj[Ye], Yv, Adj[Yv]]
    ) Yd 
  + ((4 g1^2) / 5 + 6 g2^2 - 9 trace[Yd, Adj[Yd]] - 3 trace[Ye, Adj[Ye]]) MatMul[Yd, Adj[Yd], Yd] 
  + (4 g1^2 MatMul[Yu, Adj[Yu], Yd]) / 5 
  - 3 trace[Yu, Adj[Yu]] MatMul[Yu, Adj[Yu], Yd] 
  - trace[Yv, Adj[Yv]] MatMul[Yu, Adj[Yu], Yd] 
  - 4 MatMul[Yd, Adj[Yd], Yd, Adj[Yd], Yd] 
  - 2 MatMul[Yd, Adj[Yd], Yu, Adj[Yu], Yd] 
  - 2 MatMul[Yu, Adj[Yu], Yu, Adj[Yu], Yd]
 }, 
 {
  "Yu", 
    (
    - (13 g1^2) / 15 
    - 3 g2^2 
    - (16 g3^2) / 3 
    + 3 trace[Yu, Adj[Yu]] 
    + trace[Yv, Adj[Yv]]
    ) Yu 
  + MatMul[Yd, Adj[Yd], Yu] 
  + 3 MatMul[Yu, Adj[Yu], Yu], 
    (
      (2743 g1^4) / 450 
    + g1^2 g2^2 
    + (15 g2^4) / 2 
    + (136 g1^2 g3^2) / 45 
    + 8 g2^2 g3^2 
    - (16 g3^4) / 9 
    + (4 (g1^2 + 20 g3^2) trace[Yu, Adj[Yu]]) / 5 
    - 3 trace[Yd, Adj[Yd], Yu, Adj[Yu]] 
    - trace[Ye, Adj[Ye], Yv, Adj[Yv]] 
    - 9 trace[Yu, Adj[Yu], Yu, Adj[Yu]] 
    - 3 trace[Yv, Adj[Yv], Yv, Adj[Yv]]
    ) Yu 
  + ((2 g1^2) / 5 - 3 trace[Yd, Adj[Yd]] - trace[Ye, Adj[Ye]]) MatMul[Yd, Adj[Yd], Yu] 
  + (2 g1^2 MatMul[Yu, Adj[Yu], Yu]) / 5 
  + 6 g2^2 MatMul[Yu, Adj[Yu], Yu] 
  - 9 trace[Yu, Adj[Yu]] MatMul[Yu, Adj[Yu], Yu] 
  - 3 trace[Yv, Adj[Yv]] MatMul[Yu, Adj[Yu], Yu] 
  - 2 MatMul[Yd, Adj[Yd], Yd, Adj[Yd], Yu] 
  - 2 MatMul[Yu, Adj[Yu], Yd, Adj[Yd], Yu] 
  - 4 MatMul[Yu, Adj[Yu], Yu, Adj[Yu], Yu]
 }
}
