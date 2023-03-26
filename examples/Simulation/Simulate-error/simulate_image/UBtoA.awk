#! /bin/awk -f
#
#	take a crystal orientation (U) and cell parameters (B matrix)
#	and generate a MOSFLM orientation (A) matrix 
#
#
# example input:
# CELL 74 74 34 90 90 90
# MISSET 10.0 20.0 30.0
#
# another example:
# CELL 74 74 34 90 90 90
# UMATRIX -0.84965 0.00295893 0.527339 0.479689 -0.411082 0.775184 0.219073 0.911594 0.347857
#
# random orientation example:
# CELL 74 74 34 90 90 90
# RANDOM 
#
# random orientation with "mosaic spread" of 1 degree and seed=1231321 example:
# CELL 74 74 34 90 90 90
# RANDOM 1 1231321
#
#
# multiple MISSET, RANDOM or UMATRIX enties may be made, and will be applied in turn
#
#


BEGIN{
PI = 2*atan2(1,0)
RTD = 180/PI

# default to lysozyme cell?
AX = 0.0147471  
AY = 0.000000  
AZ = 0.000000  
BX = 0.00737357
BY = 0.0127714
BZ = 0.000000
CX = 0.000000
CY = 0.000000
CZ = 0.00411184

# default rotation matrix
uxx=1;uxy=0;uxz=0; 
uyx=0;uyy=1;uyz=0;
uzx=0;uzy=0;uzz=1;

lambda = 1
}

# take pseudo-mosflm input
/^MATRIX/ && NF==10{AX=$2;BX=$3;CX=$4;AY=$5;BY=$6;CY=$7;AZ=$8;BZ=$9;CZ=$10; Amatrix=1}

# make a B-matrix out of any cell
/^CELL/ && NF==7{
    a=$2;b=$3;c=$4;alpha=$5*PI/180;beta=$6*PI/180;gamma=$7*PI/180;

    # if cell is defined, create a B-matrix
    s=(alpha+beta+gamma)/2;skew=sin(s)*sin(s-alpha)*sin(s-beta)*sin(s-gamma);
    if(skew<0) skew=-skew; if(skew==0) skew=1e-30;
    Volume = 2*a*b*c*sqrt(skew);
    a_star = b*c*sin(alpha)/Volume;
    b_star = c*a*sin(beta)/Volume;
    c_star = a*b*sin(gamma)/Volume;
    cos_alphastar = (cos(beta)*cos(gamma)-cos(alpha))/(sin(beta)*sin(gamma));
    cos_betastar  = (cos(gamma)*cos(alpha)-cos(beta))/(sin(gamma)*sin(alpha));
    cos_gammastar = (cos(alpha)*cos(beta)-cos(gamma))/(sin(alpha)*sin(beta));
    sin_alphastar = sqrt(1-cos_alphastar^2);
    sin_betastar  = sqrt(1-cos_betastar^2);
    sin_gammastar = sqrt(1-cos_gammastar^2);
    B_AX=a_star;
    B_BX=b_star*cos_gammastar;
    B_CX=c_star*cos_betastar;
    B_AY=0;
    B_BY=b_star*sin_gammastar;
    B_CY=c_star*(cos_alphastar-cos_betastar*cos_gammastar)/sin_gammastar;
    B_AZ=0;
    B_BZ=0;
    B_CZ=c_star*Volume/(a*b*c*sin_gammastar);

    matmult=1;
}

# multiply any U-matrix by the current B-matrix
/^UMAT/ && NF==10{
    uxx=$2;uxy=$3;uxz=$4; uyx=$5;uyy=$6;uyz=$7; uzx=$8;uzy=$9;uzz=$10;
    matmult=1;
}

# use misseting angles (in the same way mosflm does)
/^MISSET/ && NF>=4{
    rotX=$2/RTD; rotY=$3/RTD; rotZ=$4/RTD;
    
    # rotate around z axis
    rxx= cos(rotZ); rxy=-sin(rotZ); rxz= 0; 
    ryx= sin(rotZ); ryy= cos(rotZ); ryz= 0;
    rzx= 0;         rzy= 0;         rzz= 1;

    # multiply U = U*R
    nxx=uxx; nxy=uxy; nxz=uxz;
    nyx=uyx; nyy=uyy; nyz=uyz;
    nzx=uzx; nzy=uzy; nzz=uzz;
    uxx = nxx*rxx + nxy*ryx + nxz*rzx;
    uxy = nxx*rxy + nxy*ryy + nxz*rzy;
    uxz = nxx*rxz + nxy*ryz + nxz*rzz;
    uyx = nyx*rxx + nyy*ryx + nyz*rzx;
    uyy = nyx*rxy + nyy*ryy + nyz*rzy;
    uyz = nyx*rxz + nyy*ryz + nyz*rzz;
    uzx = nzx*rxx + nzy*ryx + nzz*rzx;
    uzy = nzx*rxy + nzy*ryy + nzz*rzy;
    uzz = nzx*rxz + nzy*ryz + nzz*rzz;

    # rotate around y axis
    rxx= cos(rotY); rxy= 0;         rxz= sin(rotY); 
    ryx= 0;         ryy= 1;         ryz= 0;
    rzx=-sin(rotY); rzy= 0;         rzz= cos(rotY);

    # multiply U = U*R
    nxx=uxx; nxy=uxy; nxz=uxz;
    nyx=uyx; nyy=uyy; nyz=uyz;
    nzx=uzx; nzy=uzy; nzz=uzz;
    uxx = nxx*rxx + nxy*ryx + nxz*rzx;
    uxy = nxx*rxy + nxy*ryy + nxz*rzy;
    uxz = nxx*rxz + nxy*ryz + nxz*rzz;
    uyx = nyx*rxx + nyy*ryx + nyz*rzx;
    uyy = nyx*rxy + nyy*ryy + nyz*rzy;
    uyz = nyx*rxz + nyy*ryz + nyz*rzz;
    uzx = nzx*rxx + nzy*ryx + nzz*rzx;
    uzy = nzx*rxy + nzy*ryy + nzz*rzy;
    uzz = nzx*rxz + nzy*ryz + nzz*rzz;

    # rotate around x axis
    rxx= 1;         rxy= 0;         rxz= 0; 
    ryx= 0;         ryy= cos(rotX); ryz=-sin(rotX);
    rzx= 0;         rzy= sin(rotX); rzz= cos(rotX);

    # multiply U = U*R
    nxx=uxx; nxy=uxy; nxz=uxz;
    nyx=uyx; nyy=uyy; nyz=uyz;
    nzx=uzx; nzy=uzy; nzz=uzz;
    uxx = nxx*rxx + nxy*ryx + nxz*rzx;
    uxy = nxx*rxy + nxy*ryy + nxz*rzy;
    uxz = nxx*rxz + nxy*ryz + nxz*rzz;
    uyx = nyx*rxx + nyy*ryx + nyz*rzx;
    uyy = nyx*rxy + nyy*ryy + nyz*rzy;
    uyz = nyx*rxz + nyy*ryz + nyz*rzz;
    uzx = nzx*rxx + nzy*ryx + nzz*rzx;
    uzy = nzx*rxy + nzy*ryy + nzz*rzy;
    uzz = nzx*rxz + nzy*ryz + nzz*rzz;

    matmult=1;
}

# apply an isotropic random orientation
/^RANDOM/{
    pi=4*atan2(1,1);

    # default "spread" to full sphere
    spread=90;
    # ... but user can specify a narrower range (in degrees)
    if($2+0>0) spread=$2;
    # default to re-seed random number generator each time, but user can specify a seed
    if($3!="") {
	srand($3);
    }else{
	srand();
    }

    # pick three random numbers on [0:1]
    r1=2*rand()-1;
    r2=2*rand()-1;
    r3=2*rand()-1;

    # use a quaternion recipie I found on the internet
    mos=spread/RTD;
    xyrad=sqrt(1-r2**2);
    rot = mos*(1-r3**2)**(1./3.);
    v1=xyrad*sin(pi*r1);v2=xyrad*cos(pi*r1);v3=r2;
    t1 =  cos(rot);
    t2 =  1 - t1;
    t3 =  v1*v1;
    t6 =  t2*v1;
    t7 =  t6*v2;
    t8 =  sin(rot);
    t9 =  t8*v3;
    t11 = t6*v3;
    t12 = t8*v2;
    t15 = v2*v2;
    t19 = t2*v2*v3;
    t20 = t8*v1;
    t24 = v3*v3;

    # now assemble rotation matrix from quaternions
    rxx = t1 + t2*t3;
    rxy = t7 - t9;
    rxz = t11 + t12;
    ryx = t7 + t9;
    ryy = t1 + t2*t15;
    ryz = t19 - t20;
    rzx = t11 - t12;
    rzy = t19 + t20;
    rzz = t1 + t2*t24;

    # multiply U = U*R
    nxx=uxx; nxy=uxy; nxz=uxz;
    nyx=uyx; nyy=uyy; nyz=uyz;
    nzx=uzx; nzy=uzy; nzz=uzz;
    uxx = nxx*rxx + nxy*ryx + nxz*rzx;
    uxy = nxx*rxy + nxy*ryy + nxz*rzy;
    uxz = nxx*rxz + nxy*ryz + nxz*rzz;
    uyx = nyx*rxx + nyy*ryx + nyz*rzx;
    uyy = nyx*rxy + nyy*ryy + nyz*rzy;
    uyz = nyx*rxz + nyy*ryz + nyz*rzz;
    uzx = nzx*rxx + nzy*ryx + nzz*rzx;
    uzy = nzx*rxy + nzy*ryy + nzz*rzy;
    uzz = nzx*rxz + nzy*ryz + nzz*rzz;

    matmult=1;
}




/^WAVE/{lambda=$2}


# multiply A = B*U if we just updated one of them
matmult==1 {
matmult=0
nxx=B_AX; nxy=B_BX; nxz=B_CX;
nyx=B_AY; nyy=B_BY; nyz=B_CY;
nzx=B_AZ; nzy=B_BZ; nzz=B_CZ;
AX = uxx*nxx + uxy*nyx + uxz*nzx;
BX = uxx*nxy + uxy*nyy + uxz*nzy;
CX = uxx*nxz + uxy*nyz + uxz*nzz;
AY = uyx*nxx + uyy*nyx + uyz*nzx;
BY = uyx*nxy + uyy*nyy + uyz*nzy;
CY = uyx*nxz + uyy*nyz + uyz*nzz;
AZ = uzx*nxx + uzy*nyx + uzz*nzx;
BZ = uzx*nxy + uzy*nyy + uzz*nzy;
CZ = uzx*nxz + uzy*nyz + uzz*nzz;
}

END {
    if(NR==0){
        print "usage:"
        print "UBtoA.awk << EOF"
	print "CELL 74 74 34 90 90 90"
	print "RANDOM"
	print "EOF"
        exit
    }
    lambda_star=1/lambda

    printf "%12.8f%12.8f%12.8f\n", AX*lambda+1e-9,BX*lambda+1e-9,CX*lambda+1e-9;
    printf "%12.8f%12.8f%12.8f\n", AY*lambda+1e-9,BY*lambda+1e-9,CY*lambda+1e-9;
    printf "%12.8f%12.8f%12.8f\n", AZ*lambda+1e-9,BZ*lambda+1e-9,CZ*lambda+1e-9;
    printf "%12.3f%12.3f%12.3f\n",0,0,0;
    printf "%12.7f%12.7f%12.7f\n", uxx,uxy,uxz;
    printf "%12.7f%12.7f%12.7f\n", uyx,uyy,uyz;
    printf "%12.7f%12.7f%12.7f\n", uzx,uzy,uzz;
    printf "%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f\n",a,b,c,alpha*RTD,beta*RTD,gamma*RTD;
    printf "%12.3f%12.3f%12.3f\n",0,0,0;

}


