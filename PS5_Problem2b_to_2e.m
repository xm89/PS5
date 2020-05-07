clear all
%% b)
ca=0:1:200;
cr1(:)=(1+100*ca.^2)./(1+ca.^2);
figure(1)
plot(ca,cr1,'k-','LineWidth',1.5)
hold on

cr2=[];
syms cr_x
for i=1:length(ca)
    e2=30*ca(i)-(100+5000*ca(i)^2)/(1+ca(i)^2+cr_x^2);
    cr_result=solve(e2,cr_x);
    check=double(cr_result);
    for j=1:length(check)
        if check(j)>0
            cr2(end+1)=real(check(j));
        end
    end
end
plot(ca(1:length(cr2)),cr2,'b-','LineWidth',1.5)
hold on
xlabel('Ca')
ylabel('Cr')
%% c)
ca=linspace(0,200,20);
cr=linspace(0,100,20);
[x_axis,y_axis]=meshgrid(ca,cr);
u=zeros(size(x_axis));
v=zeros(size(y_axis));
t=0;
for k=1:numel(x_axis)
    yprime=JSPODE(t,[x_axis(k);y_axis(k)]);
    u(k)=yprime(1);
    v(k)=yprime(2);
end

quiver(x_axis,y_axis,u,v,'r');
hold on
%% add solution for ODE
ca0=1;
cr0=10;
[t_solve,y_solve]=ode45(@JSPODE,[0,20],[ca0;cr0]);
plot(y_solve(:,1),y_solve(:,2),'r','LineWidth',1.2)
plot(y_solve(1,1),y_solve(1,2),'bo') % start point
plot(y_solve(end,1),y_solve(end,2),'ks') % end point
%% function
function dy=JSPODE(t,y)
dy(1)=-30*y(1)+(100+5000*y(1).^2)./(1+y(1).^2+y(2).^2);
dy(2)=-y(2)+(1+100*y(1).^2)./(1+y(1).^2);
dy=dy(:);
end