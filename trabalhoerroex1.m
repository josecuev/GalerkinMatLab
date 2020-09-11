%% Garlekin Method
%  HW1 - Exercicio 1

clc
clear all
close all

%% Proposta de solucao

itertations=10;
display(['Iterations: ', num2str(itertations)])
fprintf('\n%11s%11s%11s%11s%11s%11s%11s\n','N IF', 'max u', 'max du/dx', 'L2 u', 'L2 du/dx', 'Eng u', 'Eng du/dx');

for i=1:itertations
    
    
    [erronormamaxu(i), erronormaL2u(i), erronormaenergiau(i),erronormamaxdu(i), erronormaL2du(i), erronormaenergiadu(i)]=galerkinsol(i, false);
    

    x = rand(5,1);
    y = rand(5,1);
    [r,t] = cart2pol(x,y);
    fprintf('%10.0f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',[i, erronormamaxu(i),erronormamaxdu(i),erronormaL2u(i),erronormaL2du(i),erronormaenergiau(i),erronormaenergiadu(i)].');

    figure(1)
    plot(1:i,erronormamaxu,'color','#0072BD');
    hold on;
    plot(1:i,erronormaL2u,'color','#D95319');
    plot(1:i,erronormaenergiau,'color','#EDB120');
    legend('Norma max','Norma L2 ','Norma Energia' )
    title('Erro vs n�mero de fun��es de interpola��o - u(x)')

    
    figure(2);
    plot(1:i,erronormamaxdu,'color','#0072BD');
    hold on;
    plot(1:i,erronormaL2du,'color','#D95319');
    plot(1:i,erronormaenergiadu,'color','#EDB120');
    legend('Norma max','Norma L2 ','Norma Energia' )
    title('Erro vs n�mero de fun��es de interpola��o - du(x)/dx')
    drawnow;
    
end



%% Defini��o de fun��es

function [erronormamaxu,erronormaL2u, erronormaenergiau,erronormamaxdu, erronormaL2du, erronormaenergiadu]=galerkinsol(n,showplots)

    syms x g h a
    f = a*x^3;   %Fun��o f(x)
    limsup=1;    %Limite superior do intervalo

    %Fun��es Ni (e suas derivadas dNi) usadas na aproxima��o de vh
    for i=1:n
        N(i) = 1-x^i;
        dN(i) = diff(N(i),x);
    end
    
    %Fun��es Nnj (e suas derivadas dNnj) usadas na aproxima��o de gh
    Nn1 = x;
    dNn1 = diff(Nn1,x);

    %Montagem do sistema de Equa��es
    for i=1:n
        for j=1:n
            K(i,j) = int(dN(i)*dN(j),x,0,1);
        end
        F(i) = int(f*N(i),x,0,1) + subs(N(i),x,0)*h - int(dN(i)*dNn1,x,0,1)*g;
    end

    %Solu��o das constantes vi (da aproxima��o de vh):
    v = F/K;
    v=transpose(v);
    %Solu��o de vh
    vh=0;
    for i=1:n
        vh = vh + v(i)*N(i);
    end

    if showplots
        %Solu��o Exata u:
        display('************************************************')
        display('Solu��o Exata da fun��o u:')
        u = g + (1-x)*h + int( int(f,x,0,x) ,x,x,limsup); u=expand(u)
        %Solu��o Aproximada uh:
        display(' ')
        display('Solu��o Aproximada da fun��o uh:')
        uh = vh + g*Nn1;  uh=expand(uh)
        display('************************************************')
        display('Solu��o Exata da derivada da fun��o u:')
        du = diff(u,x)
        display(' ')
        display('Solu��o Aproximada da derivada da fun��o uh:')
        duh = diff(uh,x)
        display('************************************************')
    else
        %Solu��o Exata u
        u = g + (1-x)*h + int( int(f,x,0,x) ,x,x,limsup); u=expand(u);
        du = diff(u,x);
        %Solu��o Aproximada uh:
        uh = vh + g*Nn1;  uh=expand(uh);
        duh = diff(uh,x);
    end
    %% Calculo do erro

    a = 200; %Constante da fun��o f
    g = 50;  %Constante da C.C. Essencial
    h = 100; %Constante da C.C. Natural
    x=0:0.001:limsup; % Interval

    %Calculo pela norma do max
    
    Arrayu=eval(u);
    Arrayuh=eval(uh);
    variationu=abs(Arrayu-Arrayuh);
    [maxdifu,posu]=max(variationu);
    erronormamaxu=(maxdifu/abs(Arrayu(posu)));
    
    Arraydu=eval(du);
    Arrayduh=eval(duh);
    variationdu=abs(Arraydu-Arrayduh);
    [maxdifdu,posdu]=max(variationdu);
    erronormamaxdu=(maxdifdu/abs(Arraydu(posdu)));

    %Calculo pela norma L2

    syms x;
    
    erronormaL2u=double(sqrt(int(eval((u-uh)^2),x,0,limsup)/(int(eval(u^2),x,0,limsup))));
    erronormaL2du=double(sqrt(int(eval((du-duh)^2),x,0,limsup)/(int(eval(du^2),x,0,limsup))));

    %Calculo pela norma do Energia
    k=1;
    c=0;
    numerator=sqrt(0.5*int(eval(k*diff(u-uh,x)^2+c*(u-uh)^2),x,0,limsup));
    denominator=sqrt(0.5*int(eval(k*diff(u,x)^2+c*u^2),x,0,limsup));
    erronormaenergiau=double(numerator/denominator);
    
    numerator=sqrt(0.5*int(eval(k*diff(du-duh,x)^2+c*(du-duh)^2),x,0,limsup));
    denominator=sqrt(0.5*int(eval(k*diff(du,x)^2+c*du^2),x,0,limsup));
    erronormaenergiadu=double(numerator/denominator);

    %% Graficos
    % Grafico das fun��es de interpola��o

    if showplots
        x=0:0.001:limsup; % Interval
        figure
        if isempty(symvar(Nn1))==0
            plot(x,eval(Nn1),'b')
            ylim([ min([eval(Nn1),eval(N(i))]), 1.1*max([eval(Nn1),eval(N(i))]) ] )
        else
            plot(x,Nn1,'.b')
            ylim([ min([Nn1,eval(N(i))]), 1.1*max([Nn1,eval(N(i))]) ] )
        end
        hold on
        for i=1:n
            plot(x,eval(N(i)),'-r')
        end
        legend('N_{n+1}','N_i')
        title('Grafico das fun��es de interpola��o')
        
        % Grafico da solu��o exata u(x) e aproximada uh(x)
        figure
        plot(x,eval(uh),'r')
        hold on
        plot(x,eval(u),'b')
        legend('u^{h}','u')
        title('Grafico da solu��o exata u(x) e aproximada uh(x)')

        % Grafico do erro entre a solu��o exata u(x) e aproximada uh(x)
        figure
        plot(x,eval(u-uh),'b')
        legend('erro: u^{h}-u')
        title('Grafico do erro entre solu��o exata u(x) e aproximada uh(x)')

        % Grafico da primeira derivada da solu��o exata du(x) e aproximada duh(x)
        figure
        plot(x,eval(duh),'r')
        hold on
        plot(x,eval(du),'b')
        legend('du^{h}','du')
        title('Grafico das primeiras derivadas du(x) e duh(x) das solu��es')

        % Grafico do erro entre a primeira derivada da solu��o exata du(x) e aproximada duh(x)
        figure
        plot(x,eval(du-duh),'b')
        legend('erro: du^{h}-du')
        title('Grafico do erro entre as primeiras derivadas du(x) e duh(x) das solu��es')
    end

end

