breed [NICs NIC]
breed [ICs IC]

NICs-own [mode age radius PT_type change_neighbor r_p] ;;type is strictly for NIC of mode 0 (1 mutant, 2 non-muta
ICs-own [mode]
turtles-own [N_neighbors]


;; in this script have added:
;; * MADE r_p AN ATTRIBUE (important to stop error and fix it not to go beyound R_max)
;; -- pseudocode on interactions
;; * functions to compute rI and rT
;; * function to get neighboring PT and another for nieghboring IC (takes agent as input)
;; * functions to compute nPT1 and nI1 from them
;; -- in terms of simulation, same as modelNIC-copy-rayane so:
;; * only commented die to make them move (mvt is still directed, no random component)

;; ----------------------- parameters + related functions to set them -----------------------------------

globals [ p0 phi0 a b Rmax K0 K1 ;;--paramters for the simulation
         R_t delta_p delta_n nNT nNe NPT nT nI GF NF n_nonmut n_mut ] ;;--parameteres that get updated every iteration

;; --> used in setup to initialize simulation parameters (from table 3)
to initialize_parameters
  set p0 0.7      set phi0 0.06      set a 0.42
  set b 0.11      set Rmax 37.5      set K0 0.5
  set K1 0.2
end
;; --> used in go to update the time dependent params that change at each tick
to update_time_dependent_parameters
  set R_t compute_R_t     set delta_p compute_delta_p     set delta_n compute_delta_n
  set nNT countNT     set nNe countNe    set nPT countPT
  set nT nNT + nNe + nPT                 set nI countIC
  set GF nPT / nT                        set NF nNe / nT
  set n_nonmut count NICs with [mode = 1 AND PT_type = 2]
  set n_mut    count NICs with [mode = 1 AND PT_type = 1]
end

;; ------------------------------------------------------------------------------------------------------
;; ------------------------------------------------------------------------------------------------------
;; ---------------------------------- set up ------------------------------------------------------------

to setup
  clear-all
  if file-exists? "params.csv" [file-delete "params.csv"]
  initialize_parameters

  ;; --initializing all cells to NIC of mode 0
  ask patches [
    sprout-NICs 1 [
      set mode 0 set age 0 set radius distancexy 0 0 set size 0.5 set change_neighbor nobody
    ]
  ]

  ;; --setup initial PT in the center
  ask NICs [
    if radius < initial_radius [
      set mode 1
      let r random-float 1
      ifelse r < Nmm [set PT_type 2][set PT_type 1]
    ]
  ]

  ;; --setup ICs in a random corner
  setup_ICs

  color-patches-based-on-cell-type
  reset-ticks
end

;; ------------------------------------------------------------------------------------------------------
;; ------------------------------------------------------------------------------------------------------
;; --------------------------------------- go -----------------------------------------------------------

to go

  ;;-- stop simulation when all are red, no more growth
  if count NICs with [mode = 1] = 0 [ print "simulation is over"] ;;dont know how stop a forever button

  update_time_dependent_parameters ;;--rayane: this function updates variables like nNT nNe nI.. which are made global, so at each iteration they got updated automatically and can access them wherever inside go
;  let r_p 0 ;; rayane note: no need for this anymore | set as attribute

  save-params-csv "params.csv"

  ask turtles [
    let neighbor_cells turtles-on neighbors
    set N_neighbors neighbor_cells with [breed = NICs and mode = 0]
  ]

  ask NICs[
      ifelse PT_type = 1
      [ set r_p p1 radius ]
      [ set r_p p2 radius count N_neighbors ]
  ]

  ask NICs with [mode = 1 OR mode = 2][
    ifelse mode = 1 [

      ;; 1st condition: proliferate
      let r random-float 1
      let N count N_neighbors
      ifelse r_p > 0 AND r_p < r AND N > 1 [

        ;; -- chose a normal cell
        let chosen_normal_cell one-of N_neighbors

        ;; -- make 2 daughter cells
        set age 0
        let reference_PT self
        ask chosen_normal_cell [set change_neighbor reference_PT]
      ]
      [
        ifelse age > age_threshold OR radius < R_t - delta_p ;;2nd condition, -> NT
        [set mode 2]
        [set age age + 1] ;;3rd -> no change, becomes older
    ]]
    ;;else, its mode 2
    [
      if R_t - radius > delta_n + delta_p [ set mode 3 ]
      if radius > R_t - delta_p [set mode 1 let r random-float 1 ifelse r < Nmm [set PT_type 2][set PT_type 1]]
    ]
  ]


;  ;; ICs random walk
  let unbiased countT / countCells      ;;-- alpha: local tumor density (rayane - just testing)
  let biased countPT / countT           ;;-- beta:  tumor-permeable direction signal

  let maximum max list unbiased biased

  let r_walk compute_k * maximum     ;;-- commenting max ratio

  ;;-- rayane: troubleshooting
;  print (word "nT/ncell2 =" unbiased " nPT/nT =" biased)
;  print (word "r walk is :" r_walk )

  ask ICs [

    ;################## INTERACTIONS #######################
    ; pseudocode:

    ; # if IC meet PT (if there exist one PT in neighb)

    ;     -- 1: antitumor (proba rI)

    ;     # PT1 <- save the nighboring PT here (essential to save a list for 1. know to pick randomly from it and 2. count nub of successes by incrementing the number of neighb PT killed )
    ;     # randomly_chosen_PT <- choose random PT from list to move to at the end (and NK only kills one too)
    ;     # if current IC mode is 1 (CTL):
    ;     #    for each PT in PT1: set mode 4 (maybe define unstable state as mode 4)
    ;     # else (mode = 0)
    ;     #    set randomly_chosen_PT's mode to 4
    ;     # move to randomly_chosen_PT position

    ;     -- 2: neutral

    ;     # no change, stays in same state

    ;     -- 3: protumor (proba rT)
    ;     # change breed to NIC, mode to 0 (becomes empty)

    ; (considering order of 1 3 2)

    ; # else no neighboring PT
    ;    # keep on random walking

    ; ---------------------------------------------------------
    ; NOTE: can make this inside ask ICs directly, check condition if there is neighb PT (maybe when implemented compute_nPT1 check if nPT1>0)
    ;       if true -> do the interaction block
    ;       else -> do the random walk block (thats already implemented)
    ;
    ;  also note, considering implementing a function called regularize_unstable, that resets all unstable NIC states (defined maybe as state 4) to 0
    ;             this is if we wanna perform the rest in a next iteration
    ; ---------------------------------------------------------

    ;###########################################################

    let nPT1 compute_nPT1 self
    let nI1 compute_nI1 self

    ifelse nPT1 > 0 ;; --> if true perform interaction transition rules, else random walk

    [ ;; 3 interaction transition rules:

      ;; --1
      let PT1 get-PT1 self
      let random_PT_index random nPT1

    ]

    [ ;; if no neighboring PT: perform random walk
      let r random-float 1
      let chosen_normal_cell nobody
      ;;-- rayane: the problem is its always gonna be biased
      if r_walk < r
     [
        ifelse maximum = biased
        [set chosen_normal_cell min-one-of N_neighbors[radius]] ;; move toward center
        [set chosen_normal_cell one-of N_neighbors
          show True] ;; move randomly
      ]

      if chosen_normal_cell != nobody[
        let reference_IC self
        ask chosen_normal_cell [set change_neighbor reference_IC]]
    ]

  ]


  ask NICs with [mode = 0][
    if is-NIC? change_neighbor  [
      let sametype [PT_type] of change_neighbor
      set mode 1 set PT_type sametype set age 0
    ] ;; daughter cell of PT

    if is-IC? change_neighbor [
      let my-x xcor
      let my-y ycor
      let IC-x [xcor] of change_neighbor
      let IC-y [ycor] of change_neighbor
      let IC-mode [mode] of change_neighbor
      ask patch my-x my-y [sprout-ICs 1 [set mode IC-mode]]
      ask patch IC-x IC-y [
        sprout-NICs 1 [set mode 0 set radius distancexy 0 0]
        ask ICs-here [ die ]
      ]
      die
    ]
  ]
  color-patches-based-on-cell-type
  tick
end

;; ------------------------------------------------------------------------------------------------------
;; ------------------------------------------------------------------------------------------------------
;; ------------------------------------ setup ICs -------------------------------------------------------


to setup_ICs
  let ncells count patches
  let total_ICs floor (k * ncells)

  let corner random 4
  let x-min 0
  let x-max 0
  let y-min 0
  let y-max 0
  let boundary ceiling (sqrt total_ICs)

  ;; Define the bounds of the chosen corner
  if corner = 0 [ ; bottom-left
    set x-min min-pxcor
    set x-max min-pxcor + boundary
    set y-min min-pycor
    set y-max min-pycor + boundary
  ]
  if corner = 1 [ ; bottom-right
    set x-min max-pxcor - boundary
    set x-max max-pxcor
    set y-min min-pycor
    set y-max min-pycor + boundary
  ]
  if corner = 2 [ ; top-left
    set x-min min-pxcor
    set x-max min-pxcor + boundary
    set y-min max-pycor - boundary
    set y-max max-pycor
  ]
  if corner = 3 [ ; top-right
    set x-min max-pxcor - boundary
    set x-max max-pxcor
    set y-min max-pycor - boundary
    set y-max max-pycor
  ]

  ;; Get random patches from that corner
  let corner-patches patches with [
    pxcor >= x-min and pxcor <= x-max and
    pycor >= y-min and pycor <= y-max
  ]

  let chosen-patches n-of total_ICs corner-patches

  ;; Create ICs on those patches, replacing anyone there
  ask chosen-patches [
    ;; Remove any other agent on this patch
    ask turtles-here [ die ]

    ;; Create a new IC here
    sprout-ICs 1 [
      set mode random 2
    ]
  ]
end

;; ------------------------------------------------------------------------------------------------------
;; ------------------------------------------------------------------------------------------------------
;; ------------------------------------------- coloring -------------------------------------------------

to color-patches-based-on-cell-type
  ask NICs [
    set hidden? true
    if mode = 0 [ask patch-here [set pcolor 103]]
    if mode = 1 [ask patch-here [set pcolor 105]]
    if mode = 2 [ask patch-here [set pcolor 15]]
    if mode = 3 [ask patch-here [set pcolor 13]]
  ]

  ask ICs[
    set hidden? true
    ask patch-here [set pcolor yellow ]
  ]
end


;; ------------------------------------------------------------------------------------------------------
;; ------------------------------------------------------------------------------------------------------
;; --------------------------------- computing parameters -----------------------------------------------

;; ## probabilities ##
to-report p1 [r]
  let proba p0 * (1 - (r / Rmax))
  report proba
end
to-report p2 [r num_neighboring_N]
  let proba phi0 * num_neighboring_N * (1 - (r / Rmax))
  report proba
end

to-report rI_proba [a-turtle]
  let nI1 compute_nI1 a-turtle
  let nPT1 compute_nPT1 a-turtle
  let proba K0 * ( (nI1) / (nPT1) )
  report proba
end

to-report rT_proba [a-turtle] ;; --named rT_proba instead of rT bcs rT is a reserved function
  let nI1 compute_nI1 a-turtle
  let nPT1 compute_nPT1 a-turtle
  let proba K1 * ( (nPT1) / (nI1) )
  report proba
end

;; ## param on current state of the overall system ##
to-report compute_R_t
  let value 0
  let counter 0
  ask NICs with [mode = 1 ] [
    let neighbor_cells turtles-on neighbors
    let N count neighbor_cells with [breed = NICs and mode = 0]
    if N >= 1 [
      set value value + radius
      set counter counter + 1
    ]
  ]
  report value / counter
end

to-report compute_delta_n
  let value a * compute_R_t ^ (2 / 3)
  report value
end
to-report compute_delta_p
  let value b * compute_R_t ^ (2 / 3)
  report value
end

to-report compute_k
  report countIC / countCells
end


to-report countNe
  report count NICs with [mode = 3]
end
to-report countPT
  report count NICs with [mode = 1]
end
to-report countNT
  report count NICs with [mode = 2]
end
to-report countT
  report countNT + countNe + countPT
end
to-report countCells
  report count patches
end
to-report countIC
  report count ICs
end

;; ## params on one turtle ##

to-report get-PT1 [a-turtle]  ;; takes and agent and returns teh list of neighboring PT
  let neighbor_cells turtles-on [patch-here] of a-turtle  ;; gets the neighboring turtles of the given turtle
  let PT_neighbors neighbor_cells with [breed = NICs and mode = 1 ]
  report PT_neighbors
end
to-report get-I1 [a-turtle]  ;; takes and agent and returns the list of neighboring IC
  let neighbor_cells turtles-on [patch-here] of a-turtle  ;; gets the neighboring turtles of the given turtle
  let IC_neighbors neighbor_cells with [breed = ICs ]
  report IC_neighbors
end

to-report compute_nPT1 [a-turtle]
  report count get-PT1 a-turtle
end
to-report compute_nI1 [a-turtle]
  report count get-I1 a-turtle
end

;; ------------------------------------------------------------------------------------------------------
;; ------------------------------------------------------------------------------------------------------
;; ------------------------------------ saving parameters -----------------------------------------------


to save-params-csv [filename]
  ;; -- savin some params to csv
  let iteration ticks
  ifelse not file-exists? filename [
    file-open filename
    file-print "Iteration, Rt, delta p, delta n, nNT, nNe, nPT, nT, nI, nonmut, mut , GF , NF "
  ]
   [
    file-open filename
  ]

  file-print (word iteration "," R_t "," delta_p "," delta_n "," nNT "," nNe "," nPT "," nT "," nI "," n_nonmut "," n_mut ","  GF ","  NF )
  file-close
end
@#$#@#$#@
GRAPHICS-WINDOW
387
107
902
623
-1
-1
5.0
1
1
1
1
1
0
1
1
1
-50
50
-50
50
1
1
1
ticks
30.0

BUTTON
80
128
143
161
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
25
69
197
102
Nmm
Nmm
0
1
0.3
0.1
1
NIL
HORIZONTAL

INPUTBOX
93
204
159
268
initial_radius
3.0
1
0
Number

BUTTON
71
291
134
324
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
26
28
198
61
age_threshold
age_threshold
0
50
50.0
1
1
NIL
HORIZONTAL

INPUTBOX
56
341
275
401
K
0.005
1
0
Number

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.4.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
