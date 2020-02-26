# classifyForestEdge

Classify forest edges using a canopy height and canopy mask

Finds forest edges in any direction. 

Outputs raster with classified regions exposed areas [1], forested areas [2], overlapping areas [3], upwind (south-facing) [4], downwind (north-facing) [5].

Works best on 1 km tiles or 1000 x 1000 grid cells.

## inputs

'Y = matrix (2d array) of latitude (northing) in UTM coordinates'
'X = matrix (2d array) of longitude (easting) in UTM coordinates'
'YLatVec = Vector (1d array) of latitude (northing) in UTM coordinates'
'XLonVec = Vector (1d array) of longitude (easting) in UTM coordinates'

'gridSpace      = Spatial resolution of the canopy height map e.g. ASO is 3 m, NCALM is 1-m spatial resolution [m]'
'TreeHeightMultNF = How far out do you want to search as a function of tree height'
'distanceSearchSF = How far out do you want to search in m'

'CanopyMaskBox  = Canopy Mask/map (2=tree) (1=no tree) - can only have map of 1's and 2's must be same size of Y and X'
'CanopyHeight   = Canopy Height Map'
'primWindDirection = primary wind direction (degrees - use unit circle) based on a wind rose'

For example primWindDirection to 270° to get NF and SF edges. NF classified as downwind [5]. SF classified as upwind [4]. 
