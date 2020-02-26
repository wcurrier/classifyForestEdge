function [CanopyMaskOutFin] = treeHeightSingDirRegionNoOverlap_No_Upwind(Y,X,YLatVec,XLonVec,gridSpace,TreeHeightMultNF,distanceSearchSF,CanopyMaskBox,CanopyHeight,primWindDirection)
% clear all;close all;clc

% load Tuolumne_Topo_NS_Veg_NS_2.mat CanopyMask snowZ CHM
% sn=2;
% CanopyMaskBox = CanopyMask(sn).Vals;
% CanopyHeight  = CHM(sn).Vals;
% Y             = snowZ(sn).LatYMat;
% X             = snowZ(sn).LonXMat;
% YLatVec       = snowZ(sn).LatYVec;
% XLonVec       = snowZ(sn).LonXVec;
% 
% gridSpace         = 3;
% primWindDirection = 235;
% TreeHeightMult    = 3;

% Written by William Currier - currierw@uw.edu

% Find forest edges in any direction. Breaks into three regions. 
% One region is in the opposite direction of the primary direction (primWindDirection)
% Another region is in the direction of the primary direction (primWindDirection)
% Another region is areas that overlap between those two directions


% inputs

% Y = matrix of latitude (northing) in UTM coordinates
% X = matrix of longitude (easting) in UTM coordinates;

% YLatVec = Vector of latitude (northing) in UTM coordinates
% XLonVec = Vector of longitude (easting) in UTM coordinates

% gridSpace      = ASO is 3 m, NCALM is 1-m spatial resolution [m]
% distanceSearch = How far out do you want to search [m]? 3-m, 6-m, 30-m?

% CanopyMaskBox  = Canopy Mask/map (2=tree) (1=no tree) - can only have map of 1's and 2's must be same size of Y and X
% primWindDirection = primary wind direction (degrees - use unit circle) based on a wind rose

% Code

% Determine the search angles for upwind and downwind directions
upWindDirection=primWindDirection;
downWindDirection=primWindDirection+180; % flip to other side
phi=deg2rad(downWindDirection);
phi2=deg2rad(upWindDirection);

CanopyMaskBoxOut  = CanopyMaskBox;
CanopyMaskBoxOut2 = CanopyMaskBox;

                            
% Loop through the entire matrix
for r = 1:size(CanopyMaskBox,1)								% Row
	if mod(r,100) == 0
		disp(['Working on Downwind Areas row #',num2str(r)])
    end
    
        for c = 1:size(CanopyMaskBox,2)							% Column

            if CanopyMaskBox(r,c)==2

                distanceSearch=TreeHeightMultNF.*CanopyHeight(r,c); 
                if distanceSearch<4.5
                    distanceSearch=4.5;
                end
                gridCells=round(distanceSearch/gridSpace);
                

                % initialize
                xVec=nan(gridCells,1);  % potential longitude (easting) coordinates within X-m of the grid cell of interest in direction phi
                yVec=nan(gridCells,1);  % potential latitude (northing) coordinates within X-m of the grid cell of interest in direction phi

                XclosestIdx=nan(gridCells,1); % horizontal dimension indices in larger domain (before clipped)
                YclosestIdx=nan(gridCells,1); % vertical dimension indices in larger domain (before clipped)
                
                
                        % Compute easting and northing coordinates at the end of the search distance
                        xlin = distanceSearch.*cos(phi(1))+X(r,c);	% X (easting) coordinate at end of line direction phi
                        ylin = distanceSearch.*sin(phi(1))+Y(r,c);	% Y (northing) coordinate at end of line direction phi


                       if X(r,c)>=xlin && Y(r,c)<=ylin % top left quad

                            xVec(:,1)=linspace(xlin,X(r,c),gridCells).';
                            yVec(:,1)=flipud(linspace(Y(r,c),ylin,gridCells).');

                       elseif X(r,c)>=xlin && Y(r,c)>=ylin % bottom left quad

                            xVec(:,1)=linspace(xlin,X(r,c),gridCells).';
                            yVec(:,1)=linspace(ylin,Y(r,c),gridCells).';

                       elseif X(r,c)<=xlin && Y(r,c)>=ylin % bottom right quad

                            xVec(:,1)=linspace(X(r,c),xlin,gridCells).';
                            yVec(:,1)=flipud(linspace(ylin,Y(r,c),gridCells).');

                       elseif X(r,c)<=xlin && Y(r,c)<=ylin % top right quad

                            xVec(:,1)=linspace(X(r,c),xlin,gridCells).';
                            yVec(:,1)=linspace(Y(r,c),ylin,gridCells).';

                       end
                       
                       
                       
                       % Find within the larger domain the index values assosciated with the potential lat lon. Possible some are repeated depedning on how many times we iterated
                           for ii=1:gridCells
                                [XclosestIdx(ii,1),YclosestIdx(ii,1)]=findNearestIdxGridCell(YLatVec,XLonVec,yVec(ii,1),xVec(ii,1));
                           end



                       % Make a new Canopy Classification
                       YclosestIdxVec=reshape(YclosestIdx,[size(YclosestIdx,1)*size(YclosestIdx,2),1]);
                       XclosestIdxVec=reshape(XclosestIdx,[size(XclosestIdx,1)*size(XclosestIdx,2),1]);
                       YXpairs=unique([YclosestIdxVec,XclosestIdxVec],'rows'); % unique YX dimensions from the first search (4) (y first, x second) (row first, column second)

                       for nn=1:length(YXpairs(:,1))
                          if CanopyMaskBoxOut(YXpairs(nn,1),YXpairs(nn,2))==1
                                CanopyMaskBoxOut(YXpairs(nn,1),YXpairs(nn,2))=20;

                          end
                       end
            end
        end
end



% Loop through the entire matrix
for r = 1:size(CanopyMaskBox,1)								% Row
	if mod(r,100) == 0
		disp(['Working on Upwind Areas row #',num2str(r)])
    end
        for c = 1:size(CanopyMaskBox,2)							% Column

            if CanopyMaskBox(r,c)==2

%                 distanceSearch=TreeHeightMultNF.*CanopyHeight(r,c); 
                distanceSearch = distanceSearchSF;
%                 if distanceSearch<4.5
%                     distanceSearch=4.5;
%                 end
                gridCells=round(distanceSearch/gridSpace);
                
                    xVec2=nan(gridCells,1);  % potential longitude (easting) coordinates within X-m of the grid cell of interest in direction phi
                    yVec2=nan(gridCells,1);  % potential latitude (northing) coordinates within X-m of the grid cell of interest in direction phi
                    
                    XclosestIdx2=nan(gridCells,1); % horizontal dimension indices in larger domain (before clipped)
                    YclosestIdx2=nan(gridCells,1); % vertical dimension indices in larger domain (before clipped)
                    
                    
                        xlin2 = distanceSearch.*cos(phi2)+X(r,c);	% X (easting) coordinate at end of line direction phi
                        ylin2 = distanceSearch.*sin(phi2)+Y(r,c);	% Y (northing) coordinate at end of line direction phi

                       if X(r,c)>=xlin2 && Y(r,c)<=ylin2 % top left quad

                            xVec2(:,1)=linspace(xlin2,X(r,c),gridCells).';
                            yVec2(:,1)=flipud(linspace(Y(r,c),ylin2,gridCells).');

                       elseif X(r,c)>=xlin2 && Y(r,c)>=ylin2 % bottom left quad

                            xVec2(:,1)=linspace(xlin2,X(r,c),gridCells).';
                            yVec2(:,1)=linspace(ylin2,Y(r,c),gridCells).';

                       elseif X(r,c)<=xlin2 && Y(r,c)>=ylin2 % bottom right quad

                            xVec2(:,1)=linspace(X(r,c),xlin2,gridCells).';
                            yVec2(:,1)=flipud(linspace(ylin2,Y(r,c),gridCells).');

                       elseif X(r,c)<=xlin2 && Y(r,c)<=ylin2 % top right quad

                            xVec2(:,1)=linspace(X(r,c),xlin2,gridCells).';
                            yVec2(:,1)=linspace(Y(r,c),ylin2,gridCells).';

                       end
                       % Find within the larger domain the index values assosciated with the potential lat lon. Some are repeated depedning on how many times we iterated
                           for ii=1:gridCells
                                [XclosestIdx2(ii,1),YclosestIdx2(ii,1)]=findNearestIdxGridCell(YLatVec,XLonVec,yVec2(ii,1),xVec2(ii,1));
                           end                   


                       % Make a new Canopy Classification from second unit circle search
                       YclosestIdxVec2=reshape(YclosestIdx2,[size(YclosestIdx2,1)*size(YclosestIdx2,2),1]);
                       XclosestIdxVec2=reshape(XclosestIdx2,[size(XclosestIdx2,1)*size(XclosestIdx2,2),1]);
                       YXpairs2=unique([YclosestIdxVec2,XclosestIdxVec2],'rows'); % unique YX dimensions from the first search (4) (y first, x second) (row first, column second)

                       for nn=1:size(YXpairs2,1)
                               if CanopyMaskBoxOut2(YXpairs2(nn,1),YXpairs2(nn,2))==1
                                    CanopyMaskBoxOut2(YXpairs2(nn,1),YXpairs2(nn,2))=10;
                               end
                       end
            end
        end
end

CanopyMaskOutFin=CanopyMaskBoxOut.*CanopyMaskBoxOut2;
CanopyMaskOutFin(CanopyMaskOutFin==1)=1;    % Open area in both maps
CanopyMaskOutFin(CanopyMaskOutFin==4)=2;    % Forested area in both maps
CanopyMaskOutFin(CanopyMaskOutFin==10)=4;   % Upwind direction in one of the maps - open in another
CanopyMaskOutFin(CanopyMaskOutFin==20)=5;   % Downwind direction in one of the maps - open in another
CanopyMaskOutFin(CanopyMaskOutFin==200)=3;  % Downwind and Upwind in the maps
%
% smp;
% 
% figure
% pcolor(X,Y,CanopyMaskOutFin),colorbar,shading flat
% cmap=[blue;green;lightGray;red;cyan;purple];
% caxis([1,6])
% h=colorbar;
% colormap(gca,cmap)
% 
% figure
% chm2=CHM(sn).Vals;
% chm2(chm2<2)=NaN;
% pcolor(X,Y,chm2),colorbar,shading flat
% [M,N]=size(CHM(sn).Vals);
% numCells=M*N;
% CHMvec=reshape(CHM(sn).Vals,[numCells,1]);
% colormap(parula(length([2:round(max(CHMvec))-4])))
% caxis([2,round(max(CHMvec))-4])
% 
% 
% figure
% pcolor(X,Y,snowZ(sn).Vals),shading flat,colorbar
% caxis([quantile(snowZ(sn).ValsVec,0.05),quantile(snowZ(sn).ValsVec,0.95)])

% Output
% CanopyMaskOutFin
return




                   
