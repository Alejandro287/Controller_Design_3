function C1 = controlador(A,B,C)
syms s
sys = ss(A,B,C,0)
num = input('numero de eigenvalores: ');
if num==2
    eig1 = input('eigenvalor 1: ');
    eig2 = input('eigenvalor 2: ');
    ecd=det(s*eye(size(A))-A)
    ecd=sym2poly(ecd)
    ecr=(s-eig1)*(s-eig2)
    ecr=sym2poly(ecr)
    kprima=ecr-ecd
    kprima(1)=[]
    u=[B A*B]
    uinv=[1 ecd(2);0 1]
    pinv= u*uinv
    p=inv(pinv)
    k=kprima*p
    sys = ss(A-B*k,B,C,0)
    comprob=vpa(det((s*eye(size(A))-(A-B*k))))
    sys_k = ss(A-B*k,B,C,0);
    sys = ss(A,B,C,0);
    step(sys)
    hold on 
    step(sys_k)
end
if num==3
    eig1 = input('eigenvalor 1: ');
    eig2 = input('eigenvalor 2: ');
    eig3 = input('eigenvalor 3: ');
    ecd=det(s*eye(size(A))-A)
    ecd=sym2poly(ecd)
    ecr=(s-eig1)*(s-eig2)*(s-eig3)
    ecr=sym2poly(ecr)
    kprima=ecr-ecd
    kprima(1)=[]
    u=[B A*B A^2*B]
    uinv=[1 ecd(2) ecd(3);0 1 ecd(2); 0 0 1]
    pinv= u*uinv
    p=inv(pinv)
    k=kprima*p
    sys_k = ss(A-B*k,B,C,0)
    comprob=vpa(det((s*eye(size(A))-(A-B*k))))
    step(sys)
    hold on 
    step(sys_k)
end
if num == 4
   
    eig1 = input('eigenvalor 1: ');
    eig2 = input('eigenvalor 2: ');
    eig3 = input('eigenvalor 3: ');
    eig4 = input('eigenvalor 4: ');
    ecd=det(s*eye(size(A))-A)
    ecd=sym2poly(ecd)
    ecr=(s-eig1)*(s-eig2)*(s-eig3)*(s-eig4)
    ecr=sym2poly(ecr)
    kprima=ecr-ecd
    kprima(1)=[]
    u=[B A*B A^2*B A^3*B]
    uinv=[1 ecd(2) ecd(3) ecd(4);0 1 ecd(2) ecd(3); 0 0 1 ecd(2); 0 0 0 1]
    pinv= u*uinv
    p=inv(pinv)
    k=kprima*p
    sys_k = ss(A-B*k,B,C,0)
    comprob=vpa(det((s*eye(size(A))-(A-B*k))))
    step(sys)
    hold on 
    step(sys_k)
end
if num==5
    eig1 = input('eigenvalor 1: ');
    eig2 = input('eigenvalor 2: ');
    eig3 = input('eigenvalor 3: ');
    eig4 = input('eigenvalor 4: ');
    eig5 = input('eigenvalor 5: ');
    ecd=det(s*eye(size(A))-A)
    ecd=sym2poly(ecd)
    ecr=(s-eig1)*(s-eig2)*(s-eig3)*(s-eig4)*(s-eig5)
    ecr=sym2poly(ecr)
    kprima=ecr-ecd
    kprima(1)=[]
    u=[B A*B A^2*B A^3*B A^4*B]
    uinv=[1 ecd(2) ecd(3) ecd(4) ecd(5);0 1 ecd(2) ecd(3) ecd(4); 0 0 1 ecd(2) ecd(3); 0 0 0 1 ecd(2);0 0 0 0 1]
    pinv= u*uinv
    p=inv(pinv)
    k=kprima*p
    sys_k = ss(A-B*k,B,C,0)
    comprob=vpa(det((s*eye(size(A))-(A-B*k))))
    step(sys)
    hold on 
    step(sys_k)
end

end