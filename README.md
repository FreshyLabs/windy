Winds of Freshy 
------

A simplified library for interpolating winds vectors and speeds across a grid. Creates a location aware animated surface. 

### Credit 

* mostly all of this is taken from: 
  * https://github.com/cambecc/earth
  * https://github.com/cambecc/air 


1. Load data 
2. create canvas 
3. build an interpolation grid 
4. animate randomized particles on the canvas 
  - aging particles moves them according to a vector (x, y, magnitude)
5. color is defined by magnitude (the distance which each particle moves each step) 
