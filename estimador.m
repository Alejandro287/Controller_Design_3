function C1 = estimador(A,B,C)
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
    lprima=ecr-ecd
    lprima(1)=[]
    lprima=vec2mat(lprima, 2)
    lprima=lprima'
    u=[C; C*A]
    uinv=[1 0;ecd(2) 1]
    pinv= uinv*u
    p=inv(pinv)
    l=p*lprima
    comprob=vpa(det((s*eye(size(A))-(A-l*C))))
end
if num==3
    eig1 = input('eigenvalor 1: ');
    eig2 = input('eigenvalor 2: ');
    eig3 = input('eigenvalor 3: ');
    ecd=det(s*eye(size(A))-A)
    ecd=sym2poly(ecd)
    ecr=(s-eig1)*(s-eig2)*(s-eig3)
    ecr=sym2poly(ecr)
    lprima=ecr-ecd
    lprima(1)=[]
    lprima=vec2mat(lprima, 3)
    lprima=lprima'
    u=[C; C*A; C*A^2]
    uinv=[1 0 0;ecd(2) 1 0; ecd(3) ecd(2) 1]
    pinv= uinv*u
    p=inv(pinv)
    l=p*lprima
    comprob=vpa(det((s*eye(size(A))-(A-l*C))))
    
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
    lprima=ecr-ecd;
    lprima(1)=[];
    lprima=vec2mat(lprima, 4);
    lprima=lprima'
    u=[C; C*A; C*A^2; C*A^3]
    uinv=[1 0 0 0; ecd(2) 1 0 0; ecd(3) ecd(2) 1 0; ecd(4) ecd(3) ecd(2) 1]
    pinv= uinv*u
    p=inv(pinv)
    l=p*lprima
    comprob=vpa(det((s*eye(size(A))-(A-l*C))))
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
    lprim=ecr-ecd
    lprim(1)=[]
    lprima=vec2mat(lprima, 5)
    lprima=lprima'
    u=[C; C*A; C*A^2; C*A^3; C*A^4]
    uinv=[1 0 0 0 0; ecd(2) 1 0 0 0; ecd(3) ecd(2) 1 0 0; ecd(4) ecd(3) ecd(2) 1 0; ecd(5) ecd(4) ecd(3) ecd(2) 1]
    pinv= uinv*u
    p=inv(pinv)
    l=p*lprima
    comprob=vpa(det((s*eye(size(A))-(A-l*C))))
end

end