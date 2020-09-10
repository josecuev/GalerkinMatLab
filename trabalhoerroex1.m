%% Garlekin Method
%  HW1 - Exercicio 1

clc
clear all
close all

%% Proposta de solucao

iterations=9



for i=1:iterations
    
    
    [erronormamaxu(i), erronormaL2u(i), erronormaenergiau(i),erronormamaxdu(i), erronormaL2du(i), erronormaenergiadu(i)]=galerkinsol(i, false);
    
    if i==1
        figure;
    end
    plot(1:i,erronormamaxu,'-r');
    hold on;
    plot(1:i,erronormaL2u,'-b');
    plot(1:i,erronormaenergiau,'-g');
    
    plot(1:i,erronormamaxdu,'-y');
    hold on;
    plot(1:i,erronormaL2du,'-o');
    plot(1:i,erronormaenergiadu,'-p');
    legend('Norma max - u(x)','Norma L2 - u(x)','Norma Energia - u(x)','Norma max - du(x)/dx','Norma L2 - du(x)/dx','Norma Energia - du(x)/dx' )
    title('Erro vs iteracoes u(x)')
    drawnow;
    
end

N=1:iterations;

figure;
plot(N,erronormamaxu,'-r');
legend('Norma max')
title('Erro vs iteracoes - Norma max - u(x)')

figure;
plot(N,erronormaL2u,'-b');
legend('Norma L2')
title('Erro vs iteracoes - Norma L2 - u(x)')

figure;
plot(N,erronormaenergiau,'-g');
legend('Norma Energia')
title('Erro vs iteracoes - Norma Energia - u(x)')

figure;
plot(N,erronormamaxdu,'-y');
legend('Norma max')
title('Erro vs iteracoes - Norma max - du(x)/dx')

figure;
plot(N,erronormaL2du,'-o');
legend('Norma L2')
title('Erro vs iteracoes - Norma L2 - du(x)/dx')

figure;
plot(N,erronormaenergiadu,'-p');
legend('Norma Energia')
title('Erro vs iteracoes - Norma Energia - du(x)/dx')


function [erronormamaxu,erronormaL2u, erronormaenergiau,erronormamaxdu, erronormaL2du, erronormaenergiadu]=galerkinsol(n,showplots)

    syms x g h a
    f = a*x^3;   %Fun��o f(x)
    limsup=1;    %Limite superior do intervalo


    display(['Dimens�o: ', num2str(n)])

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
    erronormamaxu=(maxdifu/(Arrayu(posu)));
    
    Arraydu=eval(du);
    Arrayduh=eval(duh);
    variationdu=abs(Arraydu-Arrayduh);
    [maxdifdu,posdu]=max(variationdu);
    erronormamaxdu=(maxdifdu/(Arraydu(posdu)));

    %Calculo pela norma L2

    syms x;
    erronormaL2u=double(int(eval((u-uh)^2),x,0,limsup)/(int(eval(u^2),x,0,limsup)));
    erronormaL2du=double(int(eval((du-duh)^2),x,0,limsup)/(int(eval(du^2),x,0,limsup)));

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

