;; ========== INTERNAL GLOBAL VARIABLES =============
globals[
 R_max ;; Maximum radius of tumor
 R_t ;; Radius of tumor in step t
 W_p ;; Width of prolific crown
 R_p ;; Radius of prolific crown
 R_n ;; Radius of necrotic crown
 unit-conversion ;; unit conversion factor (mm/ud)
]

;; ========== CUSTOM PATCHES ==============
patches-own [
  r               ;; float - distance from the center to the lattice
  state           ;; int {0,...,7} - cell state
  tumor-cell?     ;; bool - if tumor-cell? then PC, QC or NC (is TC)
  age             ;; int - 0 when PC/US created, age++ with every step
  limit           ;; int - U(0,..,5) If age>limit PC->QC / US->(DC or prev-state)

;; -- Following attributes may be needed for next developments
;;  PC-signal       ;; float (0,1) - triggered in PC division. diffusable
;;  NK-signal       ;; float (0,1) - triggered in NK-PC conflict. diffusable
;;  immune-res      ;; float (0,1) - initialized in 0. Increments slowly for each TC if survives CTL-TC conflict
;;  previous-state  ;; int {0,...,7} - previous state of an US
;;  drug-res?       ;; bool - if drug-res? then DRC else DSC
]

;; ========== INITIALIZATION AND GO ============
to setup
  clear-all
  ;;; set world size
  set-patch-size 2
  let grid-half-size 100                               ;; (ud) - grid size of 200 uds
  resize-world (- grid-half-size) grid-half-size (- grid-half-size) grid-half-size

  ;;; set unit conversion factor and R_max
  let real-size 120                                    ;; (mm) - real size of the grid
  set unit-conversion real-size / (grid-half-size * 2) ;; (mm/ud)
  set R_max  38 / unit-conversion                      ;; with 200 uds of grid size, the conversion factor is ~ 0.6mm/ud (see ODD) => R_max = 38 mm ~ (38 / 0.6) uds ~ 64 uds

  ;;; set r for every cell and initialize with 1 PC and HCs
  ask patches [
    set r sqrt(pxcor * pxcor + pycor * pycor)
    create-HC ;; Every patch begins as a HC
  ]
  ask patch 0 0 [create-PC] ;; The tumor is born with 1 PC, in the center of the lattice
  reset-ticks
end

to go
  ;;; update global variables
  update-globals

  ;;; proliferative rules
  rules-PC

  ;;; quiescent rules
  rules-QC

  ;;; necrotic rules - TODO: DELETE IF USELESS
  rules-NC

  ;;; HC rules
  rules-HC


  tick
end

;;; ============= GLOBALS UPDATES ==================
to update-globals
    set R_t mean [r] of patches with [count neighbors with [state = 0] > 0 and tumor-cell?]
    set W_p a_p * R_t ^ ( 2 / 3 )
    set R_p R_t - W_p
    set R_n R_t - ((a_p + a_q) * R_t ^ ( 2 / 3 ))
end

;;; ============= 1 PROLIFERATIVE CELLS ===============

to rules-PC
  ask patches with [state = 1]
  [
    if-else divide? p_0 r
      [divide-in-W_p]                             ; Proliferate in neighbourhood W_p if possible if it can divide with probability P_div(p_0,r)
      [set age age + 1]                           ; Else its age increases
    let rpa random-init R_p (- R_p * 0.1)         ; uniform random number between R_p and R_p - 10%R_p and R_p
    if age > limit or r < rpa  [create-QC]        ; If the age has reached the limit or the cell is in QC zone, it turns to QC
  ]
end

to divide-in-W_p
  let divided false                                                 ; holder to track if PC has divided or not after the search.
  ; searches for HCs starting with those closest to it
  ; and then those farther away.
  foreach (range 1 (max list 2 W_p))                                ; I introduced max(2,W_p) to let the first PC proliferate.
  [x ->
    let W_p-neighbors neighbors in-radius x with [state = 0]        ; list of HC neighbours in radius x
    ; If exists at least 1 HC in neighbourhood
    ; it proliferates and invades HC
    if count W_p-neighbors > 0 [
      set divided true
      create-PC
      ask one-of W_p-neighbors [create-PC]
    ]
    ; If not divided age increases
    if not divided [set age age + 1]
  ]
end

; Check P_div and set if it can devide or not
to-report divide? [p radius]
  report random-float 1 <= P_div p radius
end

; Division probability fromula of a PC
to-report P_div [ p radius ]
  report p * ( 1 - radius / R_max )
end

;;; ============ 2 QUIESCENT CELLS ================
to rules-QC
  ask patches with [state = 2]
  [
    ;; MODIFICATION: QC to PC if there is a HC in its neighborhood.
    ;; avoids no proliferation and empty holes in the tumor if p_0 is too low at the start.
    if-else count neighbors with [state = 0] > 1 [create-PC]
    [;;With these changes I have achieved a more realistic effect of the quiescent crown
      ;;since randomness is introduced in the R_q and R_n parameters
      let rpa random-init R_p (R_p * 0.1)      ; uniform random number between R_p and R_p + 10%R_p
      let rna random-init R_n (- R_n * 0.1)    ; uniform random number between R_n and R_n - 10%R_n and R_n
      if r >= rpa  [create-PC]                 ; QC->PC if is is in PC zone
      if r <=  rna [create-NC]                 ; QC->NC if it is in NC zone
    ]
  ]
end
;;; ============ 3 NECROTIC CELLS - TODO: DELETE IF USELESS ================
to rules-NC
end

;;; ============ 4 NATURAL KILLER CELLS ============
to rules-NK
end

to born-NK? ;This serves to create NK cells based on the cell concentration
  let interruptor random-float 1
  if interruptor <= NK_Concentration
  [
    create-NK
  ]
end

;;; ============ 0 HEALTHY CELLS / EMPTY SPACES ================
to rules-HC
end

;;; ============== UTILS ==================

to create-HC
  set state 0
  set tumor-cell? false
  set pcolor white
end

to create-PC
  set state 1
  set tumor-cell? true
  set age 0
  set limit random-init 1 10 ; random int between 1 and 10. TODO: Talk with the team about these values
  set pcolor red
end

to create-QC
  set state 2
  set tumor-cell? true
  set pcolor pink
end

to create-NC
  set state 3
  set tumor-cell? true
  set pcolor 11
end
to create-NK
  set state 4
  set tumor-cell? false
  set pcolor blue
end


; Returns a random number between m and m + s - 1
to-report random-init [m s]
 report m + random s
end
@#$#@#$#@
GRAPHICS-WINDOW
255
13
665
424
-1
-1
2.0
1
10
1
1
1
0
1
1
1
-100
100
-100
100
0
0
1
ticks
30.0

BUTTON
11
28
74
61
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

BUTTON
77
28
140
61
NIL
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
144
28
236
61
go-forever
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

PLOT
676
15
896
209
Radius
Iterations (5 h)
Radius (mm)
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"R_t" 1.0 0 -2674135 true "" "plot R_t * unit-conversion"
"R_p" 1.0 0 -2064490 true "" "plot R_p * unit-conversion"
"R_n" 1.0 0 -7500403 true "" "plot R_n * unit-conversion"

SLIDER
10
68
182
101
p_0
p_0
0
1
0.7
0.01
1
NIL
HORIZONTAL

SLIDER
10
105
182
138
a_p
a_p
0
0.99
0.6
0.01
1
NIL
HORIZONTAL

SLIDER
10
138
182
171
a_q
a_q
0
1 - a_p
0.35
0.01
1
NIL
HORIZONTAL

MONITOR
900
15
964
60
R_t (mm)
R_t * unit-conversion
2
1
11

MONITOR
900
65
966
110
R_p (mm)
R_p * unit-conversion
2
1
11

MONITOR
900
114
966
159
R_n (mm)
R_n * unit-conversion
2
1
11

PLOT
1082
217
1282
367
N PC
Iterations (5 h)
Number of PC
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"N_PC" 1.0 0 -16777216 true "" "plot count patches with [state = 1]"

PLOT
879
217
1079
367
N QC
Iterations (5 h)
Number of QC
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"N_QC" 1.0 0 -16777216 true "" "plot count patches with [state = 2]"

PLOT
676
217
876
367
N NC
Iterations (5 h)
Number of NC
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"N_NC" 1.0 0 -16777216 true "" "plot count patches with [state = 3]"

PLOT
973
14
1222
209
Volume
Iterations (5 h)
Volume (cm^3)
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"volume" 1.0 0 -16777216 true "" "plot (4 / 3 * pi * (R_t * unit-conversion) ^ 3) / 1000"

MONITOR
900
164
967
209
V (cm^3)
(4 / 3 * pi * (R_t * unit-conversion) ^ 3) / 1000
2
1
11

SLIDER
12
176
184
209
NK_Concentration
NK_Concentration
0
0.1
0.02
0.005
1
NIL
HORIZONTAL

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
NetLogo 6.2.1
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
