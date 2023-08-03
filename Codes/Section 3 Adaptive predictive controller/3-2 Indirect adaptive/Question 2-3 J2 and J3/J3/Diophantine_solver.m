function [ G, F ] = Diophantine_solver( A , d )
F=[0 1 zeros(1,d-1)];
for i=3:length(F)
    F(i)=-A(2)*F(i-1)-A(3)*F(i-2);
end
F=F(2:end);
G=[-A(2)*F(end)-A(3)*F(end-1) -A(3)*F(end)];
end