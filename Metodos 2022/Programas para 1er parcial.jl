function biseccion(f,(a,b);tol_x=1e-10,tol_f=1e-10,nmax=100)      
        if b<a              
            a,b=b,a
        end
    Vp = zeros(nmax)
    Vf = zeros(nmax)
    Vc = zeros(nmax)
    Vr = zeros(nmax)
    fa = f(a)                    
    @assert fa*f(b)<0.     
    for i in 1:nmax
        p = a + 0.5*(b-a)
        fp = f(p)
        c = 0.5*(b-a)      ### cota error absoluto
        r = abs(c/p)       ### cota error relativo
        Vp[i] = p
        Vf[i] = fp
        Vc[i] = c
        Vr[i] = r
        if  r<tol_x && abs(fp)<tol_f   
             return (true,Vp[1:i],Vf[1:i],Vc[1:i],Vr[1:i])
        end
        if fa*fp<0.  
            b = p
        else 
                a = p
                fa = fp
        end
    end
    return false,Vp,Vf,Vc,Vr
end

***************************************************************************************************

function Newton_Raphson(f,df,p;tol_x=1e-10,tol_f=1e-10,n_max=50)
    Vp=zeros(n_max)
    Vf=zeros(n_max)
    Vc=zeros(n_max)
    Vr=zeros(n_max)
    for i in 1:n_max
        p1=p
        p+=-f(p)/df(p)         ###  p=p-f(p)/df(p)
        Vp[i]=p
        Vf[i]=f(p)
        Vc[i]=abs(p-p1)
        Vr[i]=Vc[i]/abs(p)
        if Vr[i]<tol_x && Vf[i]<tol_f         ### como entiente el programa (P ^ Q) v S ó P ^ (Q v S)
            return (true,Vp[1:i],Vf[1:i],Vc[1:i],Vr[1:i])
        end
    end
    return false,Vp,Vf,Vc,Vr
end

*****************************************************************************************************

function secante(f,(p0,p1);tol_x=1e-10,tol_f=1e-10,n_max=50)
###     if abs(f(p0))<abs(f(p1))      ### tomo a p1/ |f(p1)|>|f(p2)| para converger mas rápido?
###        p0,p1=p1,p0
###    end
    Vp=zeros(n_max)
    Vf=zeros(n_max)
    Vc=zeros(n_max)
    Vr=zeros(n_max)
    for i in 1:n_max
        fp1=f(p1)
        p2=p1-fp1*(p1-p0)/(fp1-f(p0))
        Vp[i]=p2
        Vf[i]=f(p2)
        Vc[i]=abs(p2-p1)
        Vr[i]=Vc[i]/abs(p2)
        if Vr[i]<tol_x && Vf[i]<tol_f
            return true,Vp[1:i],Vf[1:i],Vc[1:i],Vr[1:i]     ### return false,Vp,Vf,Vc,Vr     devuelve vectores de dim(n_max) y V[n]=0 para todo n>i / en la iteracion i se cumple las tolerancias
        end
        p0=p1
        p1=p2
    end
    return false,Vp,Vf,Vc,Vr
end

*******************************************************************************************************

function regula_falsi(f,(p0,p1);tol_x=1e-5,tol_f=1e-5,nmax=50)
    if p1<p0
        p0,p1=p1,p0 
    end                  
    @assert f(p0)*f(p1)<0.  
    Vp=zeros(nmax)
    Vf=zeros(nmax)
    Vc=zeros(nmax)
    Vr=zeros(nmax)
    for i in 1:nmax
        fp1=f(p1)
        p2=p1-fp1*(p1-p0)/(fp1-f(p0))
        fp2=f(p2)
        Vp[i]=p2
        Vf[i]=fp2
        Vc[i]=abs(p2-p1)
        Vr[i]=Vc[i]/abs(p2)
        if Vr[i]<tol_x && Vf[i]<tol_f
            return true,Vp[1:i],Vf[1:i],Vc[1:i],Vr[1:i]     ### return false,Vp,Vf,Vc,Vr     devuelve vectores de dim(n_max) y V[n]=0 para todo n>i / en la iteracion i se cumple las tolerancias
        end
        if fp2*fp1<0    ### la raiz está en [p1,p2]
            p0=p1
            p1=p2
        else           ### la raiz está en [p0,p2]
            p1=p2
        end
    end
    return false,Vp,Vf,Vc,Vr
end

************************************************************************************************************

function raiz_cuadratica(a,b,c)
    dis=b^2-4*a*c
     if dis>=0
        r1=(-b+sqrt(dis))/(2*a)
        r2=(-b-sqrt(dis))/(2*a)
    print("r1=$r1, r2=$r2")
    else
        d=-b/2/a 
        e=sqrt(-dis)/2/a
    print("r1=$d+i$e, r2=$d-i$e")
    end
end

**************************************************************************************************************

function tiro_oblicuo(angulo,v0;x0=0,y0=0,puntos=600)
    a=-9.8                      ### aceleracion
    v0_x=v0*cos(angulo*pi/180)        ### velodidad inicial
    v0_y=v0*sin(angulo*pi/180)
    ###r1,r2=raiz_cuadratica(a/2,v0_y,y0)
    ###return r1,r2
    dt=(r2-r1)/puntos
    t=0
    Vx=zeros(puntos)
    Vy=zeros(puntos)
    Rx=zeros(puntos)
    Ry=zeros(puntos)
    t_=zeros(puntos)
    for i in 1:puntos
    v_x=v0_x
    v_y=a*t+v0_y
    r_x=v0_x*t+x0
    r_y=a/2*t^2+v0_y*t+y0
    t=t+dt
    Vx[i]=v_x
    Vy[i]=v_y
    Rx[i]=r_x
    Ry[i]=r_y
    t_[i]=t
    end
    return Vx,Vy,Rx,Ry,t_
end

*********************************************************************

function cambio_de_base(n,b)
    @assert b>=2 && b<=16
    digitos = "0123456789ABCDEF"
    s = ""
    while n>0
    r = n%b
    s = digitos[r+1] * s      ### digito N°1--->0  
    n = n÷b
    end
    return s
end

***********************************************************************

***Comandos Plots***
using Plots
plot(a:0.01:b,f,xlabel="x",ylabel="f(x)",label="")
scatter(Vc, yscale=:log10, title="Errores vs Iteraciones", ylabel="Errores en escala logaritmica", xlabel="i")
scatter!(Vr)
scatter(0:34,Vr,xlabel="i",ylabel="Error Relatiivo",label="Biseccion",yscale=:log10,ylim=(-0.1,1.1))
p1 = plot(t,Rx,title="x(t)")
p2 = plot(t,Ry,title="y(t)")
p3 = plot(t,Vy,title="vy(t)")
p4 = plot(Rx,Ry,title="y(x)")
scatter(p1,p2,p3,p4)
# importamos el paquete Plots.
using Plots
using LaTeXStrings # Esto sirve para incorporar LaTeX al gráfico (ej. en ejes, leyendas, etc.).
# Definimos un rango de valores que comienza en -π, va de a pasos 0.01 y termina en π.
x = -π:0.01:π
# Definimos la función que queremos plotear en dicho rango
f(x) = cos(x) + 2.0*sin(x)
# Usando broadcasting (i.e. poniendo un punto detras del nombre de la función), aplicamos la función "f"
# al rango previamente definido "x" para crear un vector "y".
y = f.(x)
# Ploteamos los valores de y vs x, usando una linea punteada de color rojo, ancho 2.5 y leyenda "f".
.6plot(x,y,linestyle=:dash,linecolor=:red,linewidth=2.5,label="f")
# Agregamos al plot (por eso usamos plot! en vez de plot) título y nombres a los ejes.
plot!(xlabel=L"x",ylabel=L"y",title="Gráfico ejemplo.")
# Definimos otra curva
x = -2π:0.01:2π
g(x) = abs(cos(x))
y = g.(x)
# y la agregamos al gráfico.
plot!(x,y,linestyle=:solid,linecolor=:blue,linewidth=2.5,label="g")

*********************************************************************


