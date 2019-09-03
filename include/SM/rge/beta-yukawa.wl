{
 {
  "Yv", 
    (
    - (9 g1^2) / 20 
    - (9 g2^2) / 4 
    + 3 trace[Yd, Adj[Yd]] 
    + trace[Ye, Adj[Ye]] 
    + 3 trace[Yu, Adj[Yu]] 
    + trace[Yv, Adj[Yv]]
    ) Yv 
  - (3 (MatMul[Ye, Adj[Ye], Yv] - MatMul[Yv, Adj[Yv], Yv])) / 2, 
  (
    2 (
        21 g1^4 
      - 54 g1^2 g2^2 
      - 230 g2^4 
      + 60 lambda^2 
      + 25 (g1^2 + 9 g2^2 + 32 g3^2) trace[Yd, Adj[Yd]] 
      + 75 (g1^2 + g2^2) trace[Ye, Adj[Ye]] 
      + 85 g1^2 trace[Yu, Adj[Yu]] 
      + 225 g2^2 trace[Yu, Adj[Yu]] 
      + 800 g3^2 trace[Yu, Adj[Yu]] 
      + 15 g1^2 trace[Yv, Adj[Yv]] 
      + 75 g2^2 trace[Yv, Adj[Yv]] 
      - 270 trace[Yd, Adj[Yd], Yd, Adj[Yd]] 
      + 60 trace[Yd, Adj[Yd], Yu, Adj[Yu]] 
      - 90 trace[Ye, Adj[Ye], Ye, Adj[Ye]] 
      + 20 trace[Ye, Adj[Ye], Yv, Adj[Yv]] 
      - 270 trace[Yu, Adj[Yu], Yu, Adj[Yu]] 
      - 90 trace[Yv, Adj[Yv], Yv, Adj[Yv]]
      ) Yv 
  + (
    - 243 g1^2 
    + 45 g2^2 
    + 300 trace[Yd, Adj[Yd]] 
    + 100 trace[Ye, Adj[Ye]] 
    + 300 trace[Yu, Adj[Yu]] 
    + 100 trace[Yv, Adj[Yv]]
    ) MatMul[Ye, Adj[Ye], Yv] 
  + 279 g1^2 MatMul[Yv, Adj[Yv], Yv] 
  + 675 g2^2 MatMul[Yv, Adj[Yv], Yv] 
  - 480 lambda MatMul[Yv, Adj[Yv], Yv] 
  - 540 trace[Yd, Adj[Yd]] MatMul[Yv, Adj[Yv], Yv] 
  - 180 trace[Ye, Adj[Ye]] MatMul[Yv, Adj[Yv], Yv] 
  - 540 trace[Yu, Adj[Yu]] MatMul[Yv, Adj[Yv], Yv] 
  - 180 trace[Yv, Adj[Yv]] MatMul[Yv, Adj[Yv], Yv] 
  + 220 MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Yv] 
  - 80 MatMul[Ye, Adj[Ye], Yv, Adj[Yv], Yv] 
  - 20 MatMul[Yv, Adj[Yv], Ye, Adj[Ye], Yv] 
  + 120 MatMul[Yv, Adj[Yv], Yv, Adj[Yv], Yv]
  ) / 80
 }, 
 {
  "Yu", 
    (
    - (17 g1^2) / 20 
    - (9 g2^2) / 4 
    - 8 g3^2 
    + 3 trace[Yd, Adj[Yd]] 
    + trace[Ye, Adj[Ye]] 
    + 3 trace[Yu, Adj[Yu]] 
    + trace[Yv, Adj[Yv]]
    ) Yu 
  - (3 (MatMul[Yd, Adj[Yd], Yu] - MatMul[Yu, Adj[Yu], Yu])) / 2, 
    (
     (
       1187 g1^4 
     - 270 g1^2 g2^2 
     - 3450 g2^4 
     + 760 g1^2 g3^2 
     + 5400 g2^2 g3^2 
     - 64800 g3^4 
     + 900 lambda^2 
     + 375 (g1^2 + 9 g2^2 + 32 g3^2) trace[Yd, Adj[Yd]] 
     + 1125 (g1^2 + g2^2) trace[Ye, Adj[Ye]] 
     + 1275 g1^2 trace[Yu, Adj[Yu]] 
     + 3375 g2^2 trace[Yu, Adj[Yu]] 
     + 12000 g3^2 trace[Yu, Adj[Yu]] 
     + 225 g1^2 trace[Yv, Adj[Yv]] 
     + 1125 g2^2 trace[Yv, Adj[Yv]] 
     - 4050 trace[Yd, Adj[Yd], Yd, Adj[Yd]] 
     + 900 trace[Yd, Adj[Yd], Yu, Adj[Yu]] 
     - 1350 trace[Ye, Adj[Ye], Ye, Adj[Ye]] 
     + 300 trace[Ye, Adj[Ye], Yv, Adj[Yv]] 
     - 4050 trace[Yu, Adj[Yu], Yu, Adj[Yu]] 
     - 1350 trace[Yv, Adj[Yv], Yv, Adj[Yv]]
     ) Yu
    ) / 600 
  + (
      (
      - 43 g1^2 
      + 45 g2^2 
      - 1280 g3^2 
      + 300 trace[Yd, Adj[Yd]] 
      + 100 trace[Ye, Adj[Ye]] 
      + 300 trace[Yu, Adj[Yu]] 
      + 100 trace[Yv, Adj[Yv]]
      ) MatMul[Yd, Adj[Yd], Yu] 
    + (
        223 g1^2 
      + 675 g2^2 
      + 1280 g3^2 
      - 480 lambda 
      - 540 trace[Yd, Adj[Yd]] 
      - 180 trace[Ye, Adj[Ye]] 
      - 540 trace[Yu, Adj[Yu]] 
      - 180 trace[Yv, Adj[Yv]]
      ) MatMul[Yu, Adj[Yu], Yu] 
    + 20 (
           11 MatMul[Yd, Adj[Yd], Yd, Adj[Yd], Yu] 
         - 4 MatMul[Yd, Adj[Yd], Yu, Adj[Yu], Yu] 
         - MatMul[Yu, Adj[Yu], Yd, Adj[Yd], Yu] 
         + 6 MatMul[Yu, Adj[Yu], Yu, Adj[Yu], Yu]
         )
    ) / 80
 }, 
 {
  "Ye", 
  (
    (
    - 9 g1^2 
    - 9 g2^2 
    + 12 trace[Yd, Adj[Yd]] 
    + 4 trace[Ye, Adj[Ye]] 
    + 12 trace[Yu, Adj[Yu]] 
    + 4 trace[Yv, Adj[Yv]]
    ) Ye 
  + 6 (MatMul[Ye, Adj[Ye], Ye] - MatMul[Yv, Adj[Yv], Ye])
  ) / 4, 
  (
    2 ( 
        1371 g1^4 
      + 270 g1^2 g2^2 
      - 1150 g2^4 
      + 300 lambda^2 
      + 125 (g1^2 + 9 g2^2 + 32 g3^2) trace[Yd, Adj[Yd]] 
      + 375 (g1^2 + g2^2) trace[Ye, Adj[Ye]] 
      + 425 g1^2 trace[Yu, Adj[Yu]] 
      + 1125 g2^2 trace[Yu, Adj[Yu]] 
      + 4000 g3^2 trace[Yu, Adj[Yu]] 
      + 75 g1^2 trace[Yv, Adj[Yv]] 
      + 375 g2^2 trace[Yv, Adj[Yv]] 
      - 1350 trace[Yd, Adj[Yd], Yd, Adj[Yd]] 
      + 300 trace[Yd, Adj[Yd], Yu, Adj[Yu]] 
      - 450 trace[Ye, Adj[Ye], Ye, Adj[Ye]] 
      + 100 trace[Ye, Adj[Ye], Yv, Adj[Yv]] 
      - 1350 trace[Yu, Adj[Yu], Yu, Adj[Yu]] 
      - 450 trace[Yv, Adj[Yv], Yv, Adj[Yv]]
      ) Ye 
  + 5 (
        3 (
            129 g1^2 
          + 225 g2^2 
          - 160 lambda 
          - 180 trace[Yd, Adj[Yd]] 
          - 60 trace[Ye, Adj[Ye]] 
          - 180 trace[Yu, Adj[Yu]] 
          - 60 trace[Yv, Adj[Yv]]
          ) MatMul[Ye, Adj[Ye], Ye] 
      - 5 (
             ( 
               27 g1^2 
             - 9 g2^2 
             - 60 trace[Yd, Adj[Yd]] 
             - 20 trace[Ye, Adj[Ye]] 
             - 60 trace[Yu, Adj[Yu]] 
             - 20 trace[Yv, Adj[Yv]]
             ) MatMul[Yv, Adj[Yv], Ye] 
           + 4 (
               - 6 MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye] 
               + MatMul[Ye, Adj[Ye], Yv, Adj[Yv], Ye] 
               + 4 MatMul[Yv, Adj[Yv], Ye, Adj[Ye], Ye] 
               - 11 MatMul[Yv, Adj[Yv], Yv, Adj[Yv], Ye]
               )
          )
      )
  ) / 400
 }, 
 {
  "Yd", 
    (
    - g1^2 / 4 
    - (9 g2^2) / 4 
    - 8 g3^2 
    + 3 trace[Yd, Adj[Yd]] 
    + trace[Ye, Adj[Ye]] 
    + 3 trace[Yu, Adj[Yu]] 
    + trace[Yv, Adj[Yv]]
    ) Yd 
  + (3 (MatMul[Yd, Adj[Yd], Yd] - MatMul[Yu, Adj[Yu], Yd])) / 2, 
    (
    - (127 g1^4) / 600 
    - (27 g1^2 g2^2) / 20 
    - (23 g2^4) / 4 
    + (31 g1^2 g3^2) / 15 
    + 9 g2^2 g3^2 
    - 108 g3^4 
    + (3 lambda^2) / 2 
    + (5 (g1^2 + 9 g2^2 + 32 g3^2) trace[Yd, Adj[Yd]]) / 8 
    + (15 (g1^2 + g2^2) trace[Ye, Adj[Ye]]) / 8 
    + (17 g1^2 trace[Yu, Adj[Yu]]) / 8 
    + (45 g2^2 trace[Yu, Adj[Yu]]) / 8 
    + 20 g3^2 trace[Yu, Adj[Yu]] 
    + (3 g1^2 trace[Yv, Adj[Yv]]) / 8 
    + (15 g2^2 trace[Yv, Adj[Yv]]) / 8 
    - (27 trace[Yd, Adj[Yd], Yd, Adj[Yd]]) / 4 
    + (3 trace[Yd, Adj[Yd], Yu, Adj[Yu]]) / 2 
    - (9 trace[Ye, Adj[Ye], Ye, Adj[Ye]]) / 4 
    + trace[Ye, Adj[Ye], Yv, Adj[Yv]] / 2 
    - (27 trace[Yu, Adj[Yu], Yu, Adj[Yu]]) / 4 
    - (9 trace[Yv, Adj[Yv], Yv, Adj[Yv]]) / 4
    ) Yd 
  + (
       (
         187 g1^2 
       + 675 g2^2 
       + 1280 g3^2 
       - 480 lambda 
       - 540 trace[Yd, Adj[Yd]] 
       - 180 trace[Ye, Adj[Ye]] 
       - 540 trace[Yu, Adj[Yu]] 
       - 180 trace[Yv, Adj[Yv]]
       ) MatMul[Yd, Adj[Yd], Yd] 
     + (
       - 79 g1^2 
       + 45 g2^2 
       - 1280 g3^2 
       + 300 trace[Yd, Adj[Yd]] 
       + 100 trace[Ye, Adj[Ye]] 
       + 300 trace[Yu, Adj[Yu]] 
       + 100 trace[Yv, Adj[Yv]]
       ) MatMul[Yu, Adj[Yu], Yd] 
    + 20 (
           6 MatMul[Yd, Adj[Yd], Yd, Adj[Yd], Yd] 
         - MatMul[Yd, Adj[Yd], Yu, Adj[Yu], Yd] 
         - 4 MatMul[Yu, Adj[Yu], Yd, Adj[Yd], Yd] 
         + 11 MatMul[Yu, Adj[Yu], Yu, Adj[Yu], Yd]
         )
    ) / 80
 }
}
