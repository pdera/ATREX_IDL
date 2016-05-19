pro test_tds
 Cij=[[384.9, 100.1,  99.2,   0.0,   0.0,   0.0],$
 	  [100.1, 246.7, 104.3,   0.0,   0.0,   0.0],$
 	  [ 99.2, 104.3, 279.2,   0.0,   0.0,   0.0],$
 	  [  0.0,   0.0,   0.0,  82.1,   0.0,   0.0],$
 	  [  0.0,   0.0,   0.0,   0.0,  93.8,   0.0],$
 	  [  0.0,   0.0,   0.0,   0.0,   0.0,  96.5]]


lp=   [4.6833, 9.9226, 5.8540]

hkl0= [0.00,-2.00,1.00]
hkls=[[0.00,-2.00,1.01],[0.01,-2.00,1.00]]

print, tdscalcs(lp, cij, hkl0, hkls)

end

function tdscalcs, uc, c4r, hkl0, hkls

    ;c4r is the 4th order elasticity tensor
    ;hkl0 is Q, the Brag peak position
    ;hkls is a nx3 array, with n hkls as (Q+p)
	;uc is a 1x3 array, with a, b, c as unit cell parameters

    n = n_elements(hkls)/3 ; number of rows in hkls file

    C=c4r

    ra = 1./uc[0]
    rb = 1./uc[1]
    rc = 1./uc[2]
    ruc = [ra,rb,rc]

    Imtds=fltarr(n) ; initial value of Imtds
    p=fltarr(n,3)   ; lattice wave vector,initial value 0

    q = hkl0 # ruc
    p = (hkls ## ruc) - (hkl0 ## ruc)

    B=fltarr(3,3) ; B=Cxp, Bijk=Cijkl*pl, intermediate rank 3 tensor

    i=0 ; index
    j=0 ; index
    k=0 ; index

    m=0 ; loop number, maxium is the number of rows in hkl file

    ; Calculate Tds intensity for mth hkl

    for m=0, n-1 do $
    begin
        B=dotprod(C,p[m])
        A=fltarr(3,3)        ; Christfolle Matrix, initial value 0
        ;-------------Calculate Christoffel Matrix, A(p)=pj*Cijkl*pl
        for i=0, 2 do $
            for k=0, 2 do $
                for j=0, 2  do $
                    A[i,k]=B[i,j,k]*p[m,j]+A[i,k]

        Imtds[m]=dotprod(transpose(hkls[m]), dotprod(invert(A), hkls[m]))
        m=m+1
    endfor
    print, Imtds
    return, Imtds

    end