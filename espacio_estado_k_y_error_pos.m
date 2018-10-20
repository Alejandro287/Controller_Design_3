function C1 = integrador(A,B,C)
syms s
sys = ss(A,B,C,0)
[row,x]=size(A);
cero=zeros(row,1);
cero=vec2mat(cero,row);
cero=cero';
AP=[A cero;-C 0]
BP=[B;0]
CP=[C 0]
sys_int=ss(AP,BP,CP,0)
num = input('numero de eigenvalores: ');
if num==3
    eig1 = input('eigenvalor 1: ');
    eig2 = input('eigenvalor 2: ');
    eig3 = input('eigenvalor 3: ');
    ecd=det(s*eye(size(AP))-AP)
    ecd=sym2poly(ecd)
    ecr=(s-eig1)*(s-eig2)*(s-eig3)
    ecr=sym2poly(ecr)
    kprima=ecr-ecd
    kprima(1)=[]
    u=[BP AP*BP AP^2*BP]
    uinv=[1 ecd(2) ecd(3);0 1 ecd(2); 0 0 1]
    pinv= u*uinv
    p=inv(pinv)
    k=kprima*p
    sys_k = ss(AP-BP*k,BP,CP,0)
    comprob=vpa(det((s*eye(size(AP))-(AP-BP*k))))
    g=CP*((s*eye(3)-AP)^-1)*BP
    step(sys_k)
    
end
if num == 4
   
    eig1 = input('eigenvalor 1: ');
    eig2 = input('eigenvalor 2: ');
    eig3 = input('eigenvalor 3: ');
    eig4 = input('eigenvalor 4: ');
    ecd=det(s*eye(size(AP))-AP)
    ecd=sym2poly(ecd)
    ecr=(s-eig1)*(s-eig2)*(s-eig3)*(s-eig4)
    ecr=sym2poly(ecr)
    kprima=ecr-ecd
    kprima(1)=[]
    u=[B AP*BP AP^2*BP AP^3*BP]
    uinv=[1 ecd(2) ecd(3) ecd(4);0 1 ecd(2) ecd(3); 0 0 1 ecd(2); 0 0 0 1]
    pinv= u*uinv
    p=inv(pinv)
    k=kprima*p
     sys = ss(AP-BP*k,BP,CP,0)
    comprob=vpa(det((s*eye(size(AP))-(AP-BP*k))))
end
if num==5
    eig1 = input('eigenvalor 1: ');
    eig2 = input('eigenvalor 2: ');
    eig3 = input('eigenvalor 3: ');
    eig4 = input('eigenvalor 4: ');
    eig5 = input('eigenvalor 5: ');
    ecd=det(s*eye(size(AP))-AP)
    ecd=sym2poly(ecd)
    ecr=(s-eig1)*(s-eig2)*(s-eig3)*(s-eig4)*(s-eig5)
    ecr=sym2poly(ecr)
    kprima=ecr-ecd
    kprima(1)=[]
    u=[BP AP*BP AP^2*BP AP^3*BP AP^4*BP]
    uinv=[1 ecd(2) ecd(3) ecd(4) ecd(5);0 1 ecd(2) ecd(3) ecd(4); 0 0 1 ecd(2) ecd(3); 0 0 0 1 ecd(2);0 0 0 0 1]
    pinv= u*uinv
    p=inv(pinv)
    k=kprima*p
end

end