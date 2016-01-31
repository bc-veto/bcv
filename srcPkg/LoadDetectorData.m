function Det = LoadDetectorData(detector)
% LoadDetectorData: Return a struct containing data on the location and 
% orientation and orientation of the specified gravitational-wave detector 
% or detector site.
%
%   Det = LoadDetectorData(detector)
%
%   detector  String. Detector site of interest.  Available detectors are 
%             'ALLEGRO' (or 'AL'), 'AURIGA' (or 'AU'), 'EXPLORER' (or
%             'EX'), 'GEO' (or 'G1', 'G'), 'LHO' (or 'H1', 'H2', 'H'),
%             'LLO' (or 'L1', 'L'), 'NAUTILUS' (or 'NA'), 'NIOBE' (or
%             'NI'), 'TAMA' (or 'T1', 'T'), and 'VIRGO' (or 'V1', 'V').
%
%   Det       Struct containing data on the detector site.  The fields are:
%
%               Det.d       Antenna response matrix.  For IFOs this is
%                             Det.d = 1/2*(Det.X*Det.X'-Det.Y*Det.Y');
%                           For bars this is
%                             Det.d = Det.n*Det.n'-[1 0 0; 0 1 0; 0 0 1]/3;
%               Det.h       Detector elevation (m).
%               Det.lambda  Detector longitude (deg E of N).
%               Det.n       Bar pointing direction in Cartesian Earth-based
%                           coordinates ([] for IFOs).
%               Det.name    Name of detector (abbreviated form: one letter 
%                           for IFOs, two letters for bars)
%               Det.phi     Detector latitude (deg N).
%               Det.psi     Detector azimuth / bar pointing direction (deg 
%                           W of N, [] for IFOs).
%               Det.V       Location of vertex/beamsplitter (IFOs) or bar 
%                           in Cartesian Earth-based coordinates (m).
%               Det.X       X arm pointing direction ([] for bars).
%               Det.Y       Y arm pointing direction ([] for bars).
%               Det.Z       Right-handed direction orthogonal to the 
%                           detector X-Y plane (IFOs), or zenith sky
%                           direction at the detector (bars). 
%
% All data for interferometers follow the WGS-84 model. For bar detectors
% the data assume a spherical Earth. As a result, the antenna response
% factors (F+, Fx) for bars are only accurate to the 1% level.
%
% $Id: LoadDetectorData.m 1941 2007-11-14 05:00:57Z jrollins $


% Notes:
% ------
% For a bar detector we have an alternative (equivalent) way to compute F+,
% Fx:
%
%   h(t) = F_+ h_+ + F_x h_x 
%
% where
%
%  F_+ = sin^2(theta) * cos(2*phi)
%  F_x = sin^2(theta) * sin(2*phi)
%
%  cos(theta) = \vec{n}\cdot\vec{k}
%
% where phi is the angle between the polarization vector and the projection
% of the long axis of the detector into the transverse plane of the wave, 
% and where \vec{n} is the pointing direction of the symmetry axis of the
% cylinder and \vec{k} is the direction to to source.


%----- Guarantee that the name of the string is all upper case.
detector = upper(detector);

if (strcmp(detector,'ALLEGRO') | strcmp(detector,'AL')) 

    %----- ALLEGRO bar detector
    %
    %----- The following assume the Earth to be a sphere!
    %
    %----- Elevation h, longitude lambda, latitude phi.
    %      From Astone et al., PRD 68 022001 (2003).  
    %        Latitude 30deg27'N
    %        Longitude 268deg50'E
    %        Azimuth 40deg % W of North
    Det.h = [];
    Det.lambda = [268+50/60+0/3600];
    Det.phi = [30+27/60+0/3600];
    Det.psi = -[40+0/60+0/3600];
    %----- Spherical-Earth vertex position.
    %      Here theta, phi are location of vertex in Earth-centered
    %      coordinates (theta=0 is North pole, phi=0 is prime meridian).
    %      Here psi is the "azimuth", the angle (counterclockwise) of the 
    %      pointing vector \vec{n} from north.
    %----- Temporary variables (ordinary spherical, rad)
    theta = (90-Det.phi)*pi/180;
    phi = (Det.lambda)*pi/180;
    %----- Vector of detector zenith sky direction.
    Omega = [ cos(phi)*sin(theta); sin(phi)*sin(theta); cos(theta)];
    Det.X = [];
    Det.Y = [];
    Det.Z = Omega;
    %----- Unit vector pointing due north from the vertex, in Earth-based
    %      coordinates. 
    Vi = [-cos(phi)*cos(theta); -sin(phi)*cos(theta); sin(theta)];
    %----- Rotate counterclockwise about Omega from North by angle psi to 
    %      pointing direction \vec{n}.
    psi = Det.psi*pi/180;
    Det.n = RotateVector(Vi,Omega,psi);
    %----- Vertex ("V") position, in Cartesian Earth-based coordinates (m).
    %      Just mulitply zenith direction vector by (approximate) Earth 
    %      radius.
    Det.V = Omega*6367089;
    %----- Clean up.
    clear theta phi Omega Vi psi
    Det.name = 'AL';
    Det.d = (Det.n*Det.n' - 1/3*[1 0 0; 0 1 0; 0 0 1]) ;

elseif (strcmp(detector,'AURIGA') | strcmp(detector,'AU')) 

    %----- AURIGA bar detector
    %
    %----- The following assume the Earth to be a sphere!
    %
    %----- Elevation h, longitude lambda, latitude phi.
    %      From Laura Cadonati.  
    %        Latitude 45deg21'12"N
    %        Longitude 11deg56'54"
    %        Azimuth 44deg % EAST of North
    Det.h = [];
    Det.lambda = [11+56/60+54/3600];
    Det.phi = [45+21/60+12/3600];
    Det.psi = -[44+0/60+0/3600];
    %----- Spherical-Earth vertex position
    %      Here theta, phi are location of vertex in Earth-centered
    %      coordinates (theta=0 is North pole, phi=0 is prime meridian).
    %      Here psi is the "azimuth", the angle (counterclockwise) of the
    %      pointing vector \vec{n} from north.
    %----- Temporary variables (ordinary spherical, rad)
    theta = (90-Det.phi)*pi/180;
    phi = (Det.lambda)*pi/180;
    %----- Vector of detector zenith sky direction.
    Omega = [ cos(phi)*sin(theta); sin(phi)*sin(theta); cos(theta)];
    Det.X = [];
    Det.Y = [];
    Det.Z = Omega;
    %----- Unit vector pointing due north from the vertex, in Earth-based
    %      coordinates. 
    Vi = [-cos(phi)*cos(theta); -sin(phi)*cos(theta); sin(theta)];
    %----- Rotate counterclockwise about Omega from North by angle psi to 
    %      pointing direction \vec{n} 
    psi = Det.psi*pi/180;
    Det.n = RotateVector(Vi,Omega,psi);
    %----- Vertex ("V") position, in Cartesian Earth-based coordinates (m).
    %      Just mulitply zenith direction vector by (approximate) Earth
    %      radius. 
    Det.V = Omega*6367089;
    %----- Clean up.
    clear theta phi Omega Vi psi
    Det.name = 'AU';
    Det.d = (Det.n*Det.n' - 1/3*[1 0 0; 0 1 0; 0 0 1]) ;

elseif (strcmp(detector,'EXPLORER') | strcmp(detector,'EX')) 

    %----- EXPLORER bar detector
    %
    %----- The following assume the Earth to be a sphere!
    %
    %----- Elevation h, longitude lambda, latitude phi.
    %      From Astone et al, CQG proceedings 2006.  These supercede 
    %      earlier (incorrect) values from Coccia et al, gr-qc/0405047 
    %        Longitude 6.05deg E
    %        Latitude 46.23deg N
    %        Azimuth 39deg % EAST of North
    Det.h = [];
    Det.lambda = [6.05+0/60+0/3600];
    Det.phi = [46.23+0/60+0/3600];
    Det.psi = -[39+0/60+0/3600];
    %----- Spherical-Earth vertex position.
    %      Here theta, phi are location of vertex in Earth-centered
    %      coordinates (theta=0 is North pole, phi=0 is prime meridian).
    %      Here psi is the "azimuth", the angle (counterclockwise) of the
    %      pointing vector \vec{n} from north.
    %----- Temporary variables (ordinary spherical, rad).
    theta = (90-Det.phi)*pi/180;
    phi = (Det.lambda)*pi/180;
    %----- Vector of detector zenith sky direction.
    Omega = [ cos(phi)*sin(theta); sin(phi)*sin(theta); cos(theta)];
    Det.X = [];
    Det.Y = [];
    Det.Z = Omega;
    %----- Unit vector pointing due north from the vertex, in Earth-based
    %      coordinates. 
    Vi = [-cos(phi)*cos(theta); -sin(phi)*cos(theta); sin(theta)];
    %----- Rotate counterclockwise about Omega from North by angle psi to
    %      pointing direction \vec{n}.
    psi = Det.psi*pi/180;
    Det.n = RotateVector(Vi,Omega,psi);
    %----- Vertex ("V") position, in Cartesian Earth-based coordinates (m).
    %      Just mulitply zenith direction vector by (approximate) Earth
    %      radius. 
    Det.V = Omega*6367089;
    %----- Clean up.
    clear theta phi Omega Vi psi
    Det.name = 'EX';
    Det.d = (Det.n*Det.n' - 1/3*[1 0 0; 0 1 0; 0 0 1]) ;

elseif (strcmp(detector,'GEO') | strcmp(detector,'G1') | strcmp(detector,'G')) 

    %----- GEO 
    %
    % All values taken from 
    % http://www.geo600.uni-hannover.de/geo600/project/location.html
    %
    % See also LSD, section 9.2, "DetectorSite.h", Anderson et al 2001.
    % WARNING: In Anderson, Brady, Creighton, and Flanagan, PRD 63 042003, 
    % 2001, the x and y arms are swapped; compare ref[45] of that paper 
    % to the above web page.  (Ie, Anderson et al call the eastern mirror 
    % the y arm; for this choice \vec{x}\times\vec{y} points into the 
    % ground rather than into the sky.)  This gives an overall sign 
    % difference for GEO_d.
    %
    % ----- IFO vertex ("GEO_V") position, LSD, section 9.2.  Extraneous
    %       precision?
    % GEO_V = [  +3.85630994953e6 ;  666598.956352 ; 5.01964141692e6  ];
    % ----- Anderson et al, 2001.  Note: x,y arms switched!
    % GEO_X = [ -0.6261 ; -0.5522 ; +0.5506 ];
    % GEO_Y = [ -0.4453 ; +0.8665 ; +0.2255 ];
    %
    %------ Location of ends of each arm, and beamsplitter.
    %       Cartesian Earth-based coordinates (m)
    GEO_X_V = [ 3856044.1 ; 667117.7 ; 5019776.0 ];  %-- eastern mirror
    GEO_Y_V = [ 3855935.6 ; 666266.5 ; 5019971.0 ];  %-- northern mirror
    GEO.V =   [ 3856311.2 ; 666597.8 ; 5019640.6 ];  %-- beam splitter
    %----- Arm pointing directions, and zenith vector.  
    %      Note {X}cross{Y} is not a unit vector, since {X}dot{Y}!=0.
    GEO_X = (GEO_X_V-GEO.V);
    GEO.X = GEO_X/(sum(GEO_X'*GEO_X))^0.5;
    GEO_Y = (GEO_Y_V-GEO.V);
    GEO.Y = GEO_Y/(sum(GEO_Y'*GEO_Y))^0.5;
    GEO_Z = cross(GEO_X,GEO_Y,1);
    GEO.Z = GEO_Z/(sum(GEO_Z'*GEO_Z))^0.5;
    %
    GEO.d = 0.5*(GEO.X*GEO.X' - GEO.Y*GEO.Y') ;
    GEO.n = [];
    GEO.psi = [];
    %----------- WGS-84 geodetic coordinates (m/deg/deg) of beam splitter
    %            Elevation h, longitude lambda, latitude phi.
    %            From uni-hannover web page.
    GEO.h      = 114.425;
    GEO.lambda = [ 9 , 48 , 25.894]*[1 ; 1/60; 1/3600];
    GEO.phi    = [52 , 14 , 42.528]*[1 ; 1/60; 1/3600];
    %----- Clean up.
    clear GEO_X_V GEO_Y_V GEO_X GEO_Y GEO_Z
    GEO.name = 'G';
    Det = GEO;

elseif (strcmp(detector,'LHO') | strcmp(detector,'H1') | ... 
        strcmp(detector,'H2') | strcmp(detector,'H')) 

    %----- LHO
    %
    % Values from Althouse et al, LIGO-P000006-D-E
    % These give d-matrix values different from those from Allen 
    % in second or third decimal place.
    LHO.X = [ -0.223891216 ; 0.799830697 ; 0.556905359 ];
    LHO.Y = [ -0.913978490 ; 0.0260953206 ; -0.404922650 ];
    LHO.Z = [ -0.338402190 ; -0.599658144 ; 0.725185541 ];
    %
    LHO.d = 0.5*(LHO.X*LHO.X' - LHO.Y*LHO.Y') ;
    %----- IFO vertex ("V") position
    %----------- Cartesian Earth-based coordinates (m)
    LHO.V = [  -2.161414928e6 ; -3.834695183e6 ; 4.600350224e6  ];
    %----------- WGS-84 geodetic coordinates (m/deg/deg)
    %            Elevation h, longitude lambda, latitude phi
    LHO.h = 142.555;
    LHO.lambda = -[119+24/60+27.565681/3600];
    LHO.phi = 46+27/60+18.527841/3600;
    LHO.fs = 16384; % Hz
    LHO.name = 'H';
    LHO.n = [];
    LHO.psi = [];
    Det = LHO;

elseif (strcmp(detector,'LLO') | strcmp(detector,'L1') | strcmp(detector,'L')) 

    %----- LLO 
    %
    % Values from Althouse et al, Review of Scientific Instruments 72, 3086 2001,
    % LIGO-P000006-D-E.
    % These give d-matrix values different from those from Allen 
    % in third decimal place.
    LLO.X = [ -0.954574615 ; -0.141579994 ; -0.262187738 ];
    LLO.Y = [ 0.297740169 ; -0.487910627 ; -0.820544948 ];
    LLO.Z = [ -0.011751435 ; -0.861335199 ; 0.507901150 ];
    %
    LLO.d = 0.5*(LLO.X*LLO.X' - LLO.Y*LLO.Y') ;
    %----- IFO vertex ("V") position
    %----------- Cartesian Earth-based coordinates (m)
    LLO.V = [  -7.427604192e4 ;  -5.496283721e6 ; 3.224257016e6  ];
    %----------- WGS-84 geodetic coordinates (m/deg/deg)
    %            Elevation h, longitude lambda, latitude phi
    LLO.h = -6.574;
    LLO.lambda = -[90+46/60+27.265294/3600];
    LLO.phi = 30+33/60+46.419531/3600;
    LLO.fs = 16384; % Hz
    LLO.name = 'L';
    LLO.n = [];
    LLO.psi = [];
    Det = LLO;

elseif (strcmp(detector,'NAUTILUS') | strcmp(detector,'NA')) 

    %----- NAUTILUS bar detector

    %----- Elevation h, longitude lambda, latitude phi.
    %      From Astone et al, gr-qc/0210053 
    %        Longitude 12.67 deg E
    %        Latitude  41.82 deg N
    %        Azimuth   44    deg EAST of North
    Det.lambda = [12.67+0/60+0/3600];
    Det.phi = [41.82+0/60+0/3600];
    Det.psi = -[44+0/60+0/3600];

    %----- The following assume the Earth to be a sphere!
    %
    %----- Spherical-Earth vertex position
    %      Here theta, phi are location of vertex in Earth-centered coordinates
    %      (theta=0 is North pole, phi=0 is prime meridian).
    %      Here psi is the "azimuth", the angle (counterclockwise) of the pointing 
    %      vector \vec{n} from north.
    %----- Temporary variables (ordinary spherical, rad)
    theta = (90-Det.phi)*pi/180;
    phi = (Det.lambda)*pi/180;
    %----- Vector of detector zenith sky direction.
    Omega = [ cos(phi)*sin(theta); sin(phi)*sin(theta); cos(theta)];
    Det.X = [];
    Det.Y = [];
    Det.Z = Omega;
    %----- Unit vector pointing due north from the vertex, in Earth-based coordinates. 
    Vi = [-cos(phi)*cos(theta); -sin(phi)*cos(theta); sin(theta)];
    %----- Rotate counterclockwise about Omega from North by angle psi to pointing 
    %      direction \vec{n} 
    psi = Det.psi*pi/180;
    Det.n = RotateVector(Vi,Omega,psi);
    %----- Vertex ("V") position, in Cartesian Earth-based coordinates (m).
    %      Just multiply zenith direction vector by (approximate) Earth radius. 
    Det.V = Omega*6367089;
    Det.name = 'NA';
    Det.d = (Det.n*Det.n' - 1/3*[1 0 0; 0 1 0; 0 0 1]) ;
    Det.h = [];

elseif (strcmp(detector,'NIOBE') | strcmp(detector,'NI')) 

    %----- NIOBE bar detector
    %
    %----- The following assume the Earth to be a sphere!
    %
    %----- Elevation h, longitude lambda, latitude phi.
    %      From Astone et al., PRD 68 022001 (2003).  
    %        Latitude -31deg56'N
    %        Longitude 115deg49'E
    %        Azimuth 0deg % EAST of North
    Det.h = [];
    Det.lambda = [115+49/60+0/3600];
    Det.phi = -[31+56/60+0/3600];
    Det.psi = [0+0/60+0/3600];
    %----- Spherical-Earth vertex position
    %      Here theta, phi are location of vertex in Earth-centered coordinates
    %      (theta=0 is North pole, phi=0 is prime meridian).
    %      Here psi is the "azimuth", the angle (counterclockwise) of the pointing 
    %      vector \vec{n} from north.
    %----- Temporary variables (ordinary spherical, rad)
    theta = (90-Det.phi)*pi/180;
    phi = (Det.lambda)*pi/180;
    %----- Vector of detector zenith sky direction.
    Omega = [ cos(phi)*sin(theta); sin(phi)*sin(theta); cos(theta)];
    Det.X = [];
    Det.Y = [];
    Det.Z = Omega;
    %----- Unit vector pointing due north from the vertex, in Earth-based coordinates. 
    Vi = [-cos(phi)*cos(theta); -sin(phi)*cos(theta); sin(theta)];
    %----- Rotate counterclockwise about Omega from North by angle psi to pointing 
    %      direction \vec{n} 
    psi = Det.psi*pi/180;
    Det.n = RotateVector(Vi,Omega,psi);
    %----- Vertex ("V") position, in Cartesian Earth-based coordinates (m).
    %      Just mulitply zenith direction vector by (approximate) Earth radius. 
    Det.V = Omega*6367089;
    %----- Clean up
    clear theta phi Omega Vi psi
    Det.name = 'NI';
    Det.d = (Det.n*Det.n' - 1/3*[1 0 0; 0 1 0; 0 0 1]) ;

elseif (strcmp(detector,'TAMA') | strcmp(detector,'T1') | strcmp(detector,'T')) 

    %----- TAMA
    % TAMA_h, TAMA_lambda, TAMA_phi are taken from Anderson, Brady, Creighton, 
    % and Flanagan, PRD 63 042003, 2001, who  
    % in turn got them via private communication with M.-K. Fujimoto.
    % These should be more accurate than the numbers in Allen, gr-qc/9607075.
    %----------- WGS-84 geodetic coordinates (m/deg/deg) of vertex.
    %            Elevation h, longitude lambda, latitude phi
    Det.h = 90;  
    Det.lambda = [139 32 9.8]*[1; 1/60; 1/3600];
    Det.phi = [35 40 35.6]*[1; 1/60; 1/3600];
    %
    %----- From Anderson et al 2001, who calculated them using the above 
    %      data and the WGS-84 model.
    %      Since the data are only given to 4 significant digits, 
    %      remormalize the units vectors so that d_ab will be traceless.
    Det.V = 1e6*[-3.946409 ; 3.366259 ; 3.699151];
    Det.X = [ 0.6490 ; 0.7608 ; 0.0];
    Det.X = Det.X/(Det.X'*Det.X)^0.5;
    Det.Y = [ -0.4437 ; 0.3785 ; -0.8123];
    Det.Y = Det.Y/(Det.Y'*Det.Y)^0.5;
    Det.Z = cross(Det.X,Det.Y);
    %----- Response matrix.
    Det.d = 0.5*(Det.X*Det.X' - Det.Y*Det.Y') ;
    Det.name = 'T';
    Det.n = [];
    Det.psi = [];
    %
    % The following assume the Earth to be a sphere!
    % Here theta, phi are location of vertex in Earth-centered coordinates
    % (theta=0 is North pole, phi=0 is prime meridian).
    % Here psi is the "azimuth", the angle (counterclockwise) of X arm from north.
    %----- Temporary variables (ordinary spherical, rad)
    %theta = (90-TAMA_phi)*pi/180;
    %phi = (TAMA_lambda)*pi/180;
    %----- Vector of detector zenith sky direction.
    %Omega = [ cos(phi)*sin(theta); sin(phi)*sin(theta); cos(theta)];
    %----- Unit vector pointing due north from the vertex, in Earth-based coordinates. 
    %Vi = [-cos(phi)*cos(theta); -sin(phi)*cos(theta); sin(theta)];
    %----- Rotate counterclockwise about Omega from North to direction of x arm 
    %      (angle psi -- just happens to be exactly 90 degrees!)
    %psi = 90*pi/180;
    %TAMA_X = RotateVector(Vi,Omega,psi);
    %----- Rotate counterclockwise about Omega by 180 deg to direction of y arm 
    %psi = 180*pi/180;
    %TAMA_Y = RotateVector(Vi,Omega,psi);
    %----- IFO vertex ("V") position
    %TAMA_V = [  -3946415.50943445 ; 3365795.25310062 ; 3699410.05882902  ];
    %Det = TAMA;

elseif (strcmp(detector,'VIRGO') | strcmp(detector,'V1') | strcmp(detector,'V')) 

    %----- Virgo
    %
    % Data are taken from Anderson, Brady, Creighton, 
    % and Flanagan, PRD 63 042003, 2001.
    % Since the data are only given to 4 significant digits, 
    % remormalize the units vectors so that d_ab will be traceless.
    Det.X = [ -0.7005 ;  0.2085 ;  0.6826 ]; 
    Det.X = Det.X/(Det.X'*Det.X)^0.5;
    Det.Y = [ -0.0538 ; -0.9691 ;  0.2408 ];
    Det.Y = Det.Y/(Det.Y'*Det.Y)^0.5;
    Det.Z = cross(Det.X,Det.Y);
    %
    Det.d = 0.5*(Det.X*Det.X' - Det.Y*Det.Y') ;
    %
    %----- IFO vertex ("V") position
    %----------- Cartesian Earth-based coordinates (m)
    Det.V = 1e6*[ 4.546374 ; 0.842990 ; 4.378577 ];
    %----------- WGS-84 geodetic coordinates (m/deg/deg)
    %            Elevation h, longitude lambda, latitude phi
    Det.h = 51.884;
    Det.lambda = [10+30/60+16.1878/3600];
    Det.phi = [43+37/60+53.0921/3600];
    % psi_x =  70.5674  % N of E  arm pointing directions
    % psi_y = 160.5674  % N of E
    Det.name = 'V';
    Det.n = [];
    Det.psi = [];

else
    
    error(['Detector ' detector ' not recognised.'])
    
end

%----- Done
return
