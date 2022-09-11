
 var res = document.getElementById("resultado");

 // ALGEBRA LINEAR -----------------------------------------------------------------------------------------------------------------------------------------------

 // PROGRESSOES ARITMETICAS E GEOMETRICAS (TERMO GERAL E SOMA DOS TERMOS) ----------------------------------------------------------------------------------------------------------------------------------------------------

 
 function enesimoTermoPA(ak,n,r,k){

  return ak+(n-k)*r;

  /*
  ak = kesimo termo
  n = posicao do termo a ser descoberto
  r = razao
  k = posicao do kesimo termo
  */
 }
 
 function somaPA (a1,an,n){
  return ((a1+an)*n)/2;

  /*
  a1 = primeiro termo
  an = enesimo termo
  n = numero de termos a serem somados
  */
 }
 
 function enesimoTermoPG(ak,n,q,k){
   return ak*Math.pow(q,(n-k));

   /*
   ak = kesimo termo
   q = razao
   n = posicao do termo a ser descoberto
   k = posicao do kesimo termo
   */
 }

 function somaPGfinita (a1,q,n){
  return (a1*(Math.pow(q,n)-1)) / (q-1);

  /*
  a1 = primeiro termo
  q = razao
  n = numero de termos a serem somados
  */

 }
 
 function somaPGinfinita(a1,q){
   return a1/(1-q);

   /*
   a1 = primeiro termo
   q = razao
   */
 }

 // GEOMETRIA PLANA (AREAS DO TRIANGULO, RETANGULO, TRAPEZIO, QUADRADO, LOSANGO, HEXAGONO E CIRCULO, PERIMETROS E APOTEMAS) --------------------------------------------------------------------------------

 function areaTriangulo(b,h){
  return (b*h)/2;
}

 
 function areaTriangEqui(l){
    return (Math.pow(3,(1/2)) / 4) * Math.pow(l,2);
 }


 function areaRetangulo(b,h){
   return (b*h);

   /*
   b = base
   h = altura
   */
 }


 function areaTrapezio(B,b,h){
   return ((B+b)*h)/2;

   /*
   B = base maior
   b = base menor
   h = altura
   */
 }


 function areaQuadrado(l){

   return (Math.pow(l,2));
 }


 function areaLosango(D,d){
   return (D*d)/2;

   /*
   D = diagonal maior
   d = diagonal menor
   */
 }


function areaPentagono(l){
    return(1/4*6.88190960236*(Math.pow(l,2)));
}


 function areaHexagono(l){
  return (3*(Math.pow(l,2)) * (Math.pow(3,(1/2)))) / 2;
 }

 
 function areaCirculo(r){
   return Math.PI*(Math.pow(r,2));
 }
 

 function perimetroCirculo(r){
   return 2*Math.PI*r;
 }
 
 // r = raio


 function perimetroGeral(n, l){
   return n*l;

   /*

   n = numero de lados
   l = valor do lado

   */
 }

 function areaGeral(p,a){
  return (p*a)/2;

  /*
  p = perimetro
  a = apotema

  */
}

function apotemaTriangulo(h){
  return (h/3);

  // h = altura
}

function apotemaQuadrado(l){
  return (l/2);
  // l = lado
}

function apotemaHexagono(l){
  return (l * Math.pow(3,(1/2)))/2;
  // l = lado
}
 
 // GEOMETRIA ESPACIAL (FORMULA DE EULER, VOLUME E AREA DA SUPERFICIE DOS PRISMAS, PIRAMIDES, CILINDROS, ESFERA E CONE (GERATRIZ)) --------------------------------------------------------------------------------


 function faceEuler(a,v){
    return (a+2-v);
 }
 
 function arestaEuler(v,f){
    return (v+f-2);
 }

 function verticeEuler(a,f){
    return (a+2-f);
 }

 /*
 f = faces
 v = vertices
 a = arestas
 */


 function volumePrisma(Ab, H){
    return (Ab*H);

    /*
    Ab = area da base
    H = altura
    */
 }

 function areaLateral(p, H){
   return (p*H);
   /*
   p = perimetro da base
   H = altura
   */
 }

function areaPrisma(Ab, Al){
  return (2*Ab+Al)
  /*
  Ab = area da base
  Al = area lateral
  */
}

function areaCuboide(a,b,c){
  return (2*(a*b+b*c+a*c));
}

function volumeCuboide(a,b,c){
  return (a*b*c);
}

function diagonalCuboide(a,b,c){
  return (Math.sqrt(Math.pow(a,2)+Math.pow(b,2)+Math.pow(c,2)));
}

 /*
  a = lado 1
  b = lado 2
  c = lado 3
  */

function alturaTetraedro(l){
  return (l*Math.sqrt(6))/3;
}

function volumeTetraedro(Ab,l){
  return (1/3*Ab*alturaTetraedro(l));
  /*
  Ab = area da base
  H = altura
  */
}

function volumePiramide(Ab,H){
  return(1/3*Ab*H);
  /*
  Ab = area da base
  H = altura
  */
}

function geratrizCone(h,r){
  return (Math.sqrt(Math.pow(h,2)*Math.pow(r,2)));
}

function areaCone(h,r){
  return Math.PI*r*geratrizCone(h,r);
}

/*
h = altura
r = raio
*/

function volumeEsfera(r){
  return (4/3*Math.PI*Math.pow(r,3));
  // r = raio
}

function areaEsfera (r){
  return (4*Math.PI*Math.pow(r,2));
  // r = raio
}

 // ESTATISTICA (ROL, MEDIA, MEDIANA, VARIANCIA E DESVIO PADRAO) ----------------------------------------------------------------------------------------------------------------------------------------------------------------
 
 function estatistica(){
        
  var salvar = parseInt(0);
  var populacao = new Array();
  var resultado = new Array();
  var limite = prompt("Digite o tamanho da população (quantidade de números)");
  
  for(var i = 0; i < limite; i++){
     var xi = parseInt((prompt("Digite um elemento")));
     populacao.push(xi);
     var salvar = salvar + xi;

  }
  
  populacao.sort(function(a, b){return a-b});
 
  var media = (salvar/limite);
  var salvar = parseInt(0);
  

  populacao.forEach(element => resultado.push(Math.pow((element-media),2))
  );

  resultado.forEach(element => salvar = salvar + element
  );


  if(limite % 2 != 0){
    var posicao = Math.ceil((limite/2));
    var mediana = populacao[posicao-1];
  }else{
    var posicao = (limite/2);
    var mediana = (populacao[posicao-1]+populacao[posicao])/2;
  }
}

/*
Rol = populacao
Media = media
Mediana = mediana
Variancia = salvar/limite
Desvio padrao = sqrt(variancia)
*/

 // ANALISE COMBINATORIA (FATORIAL, COMBINACAO, ARRANJO E PERMUTACAO) ----------------------------------------------------------------------------------------------------------------------------------------------------------------
 
 function fat(num) {
    
      for (var i = num - 1; i >= 1; i--) {
      num *= i;
    }
    return(num);
  }

  
function combinacao(n,p) {
    return (fat(n)) / (fat(n-p)*fat(p));
    /*
    n = numero de elementos
    p = tamanho dos subconjuntos
    */
}


function arranjoSimples(n,p){
    return (fat(n)) / (fat(n-p));

    /*
    n = numero de elementos
    p = tamanho dos subconjuntos
    */
}


function permutacaoSimples(n){
   return fat(n);

   /*
   n = numero de elementos
   */
}

// PROBABILIDADE (PROBABILIDADE GERAL, INTERSECCAO, UNIAO E PROBABILIDADE BINOMIAL) ----------------------------------------------------------------------------------------------------------------------------------------------------------------

function probabilidade(e,s){
  return (e/s);

  /*
  e = eventos favoraveis
  s = tamanho do espaco amostral
  */
}

function probInter(Pa,Pb){
  return (Pa*Pb);
}

function probUniao(Pa,Pb){
 
  return (Pa+Pb) - probInter(Pa,Pb);
}

 /*
  Pa = probabilidade do evento A
  Pb = probabilidade do evento B
  */


function probBinomial(n,k,p){
     return (combinacao(n,k)) * (Math.pow(p,k)) * (Math.pow((1-p),(n-k)));

/*
n = numero de vezes que o evento se repete
k = numero de sucessos
p = probabilidade de sucesso
*/

}


// TRIGONOMETRIA (SENO, COSSENO, TANGENTE, LEI DOS SENOS, LEI DOS COSSENOS, ANGULO ENTRE PONTEIROS DE RELOGIO E CONVERSAO DE RADIANOS/GRAUS) ------------------------------------------------------------------------


 function sin(x){
  return Math.sin(x*(Math.PI/180));
}

function cos(x){
  return Math.cos(x*(Math.PI/180));
}

function tan(x){
  return Math.tan(x*(Math.PI/180));
}


function leiSenos(b, sinB, sinA){
  return (b*sinA)/(sinB);
}

function leiCossenos(b,c,alfa){
  return (Math.sqrt(b*b+c*c-2*b*c*(cos(alfa))));
}

function anguloMenorPonteiros(h, min){
  
  let x = Math.abs((60*h-11*min))/2;

  if (x > 180){
    return (360-x);
  }else{
    return x;
  }
}

function tipoGrau(x){
  if (x = 90){
    return "reto";
  }else if (x = 180){
    return "raso";
  }else if (x>90){
  return "obtuso";
  }else if (x<90){
  return "agudo";
  }
}

function grausParaRadianos(alfa){
  return alfa*(Math.PI/180);
}

function radianosParaGraus(alfa){
  return alfa*(180/Math.PI);
}

function comprimentoArco(alfa, r){
  return alfa*r*(Math.PI/180);
}

// GEOMETRIA ANALITICA DE PONTO E RETA (TEOREMA DE PITAGORAS, DISTANCIA PONTO PONTO, PONTO MEDIO, DISTANCIA PONTO RETA, ALINHNAMENTO DE 3 PONTOS, EQUACAO DA RETA DADO 2 PONTOS OU 1 PONTO E 1 ANGULO, RETAS PARALELAS/PERPENDICULARES, BARICENTRO DE TRIANGULO) --------------------------------------------------------------------------


 function Pitagoras(b,c){
  return Math.sqrt(Math.pow(b,2)+Math.pow(c,2));
}
 
 function distPontoPonto(x1,x2,y1,y2){

  let distancia = Math.sqrt(Math.pow((x2-x1),2)+Math.pow((y2-y1),2));
  let pontomedio = ((x1+x2)/2)+", "+(y1+y2)/2;
  let declive = (y2-y1)/(x2-x1);
  let angulo = Math.atan(declive)*(180/Math.PI);

  return ("Distância: "+distancia+"<br> Ponto Médio: "+pontomedio+"<br> Declive: "+declive+"<br> Ângulo: "+angulo+"°");
 }


 function distPontoReta(a,b,c,x0,y0){

  let numerador = Math.abs(a*x0+b*y0+c);
  let denominador = Pitagoras(a,b);

  return numerador/denominador;
 }


 function alinhamentoPontos(xa,ya,xb,yb,xc,yc){

  let det = (xa*yb+ya*xc+xb*yc-xc*yb-yc*xa-xb*ya);
  if (det == 0){
    return "Os pontos não formam um triângulo";
  }else{
    return "Os pontos formam um triângulo. Sua área é de: "+(1/2*Math.abs(det));
  }
  }

 function equacaoReta2Pontos(x1,y1,x2,y2){
 
  let a = y1-y2;
  let b = x2-x1;
  let c = x1*y2-x2*y1;
  return a+"x + ("+b+"y) + ("+c+")";
 }


 function equacaoReta1Ponto1Angulo(xa,ya,a){
   return "y = "+tan(a)+"x + ("+(tan(a)*xa*-1+ya)+")";
 }

 function retaParalelaReduzida(m,x0,y0){
  
  let coefangular=(y0-m*x0);
  return "y = "+m+"x + ("+coefangular+")";
 }

 function retaParalelaGeral(a,b,x0,y0){

  let coefangular = (-1*(a*x0+b*y0));
  return a+"x + ("+b+"y)"+" + ("+coefangular+")"+" = 0";  
 }

 function retaPerpReduzida(m,x0,y0){
  let coefangular=(1/m*x0+y0);
  return "y = "+(-1/m)+"x + ("+coefangular+")";

 }

 function retaPerpGeral(a,b,x0,y0){
 
  let coefangular = (-1*(b*x0-a*y0));
  return b+"x + ("+(a*-1)+"y)"+" + ("+coefangular+")"+" = 0";
 }


function baricentroTriangulo(xa,ya,xb,yb,xc,yc){
  return (xa+xb+xc/3)+", "+(ya+yb+yc/3);
}

// GEOMETRIA ANALITICA DE CIRCUNFERENCIA (EQUACAO DA CIRCUNFERENCIA, POSICAO RELATIVA ENTRE PONTO E CIRCUNFERENCIA, RETA E CIRCUNFERENCIA E 2 CIRCUNFERENCIAS)


function equacaoCircunferencia(xc,yc,d){
  
  let centro = (-1*xc)+", "+(-1*yc);
  let raio = Math.sqrt(d);

  return ("Centro: "+centro+"<br> Raio: "+raio);
}

function posRelativaPontoCirc(xc,yc,d,x0,y0){
 
  let calculo = (Math.pow((x0+xc),2)+Math.pow((y0+yc),2));
  if (calculo > d){
    return "O ponto é externo";
  }else if(calculo < d){
    return "O ponto é interno";
  }else if(calculo == d){
    return "O ponto pertence a circunferência";
  }
}

function posRelativaRetaCirc(a,b,c,x0,y0,r){


  let ds = distPontoReta(a,b,c,x0,y0);
  if (ds > r){
    return "A reta não toca a circunferência";
  }else if(ds < r){
    return "A reta toca a circunferência em 2 pontos";
  }else if(ds == r){
    return "A reta toca a circunferência em 1 ponto";
  }
  
}

function posRelativaCircCirc(xc1,yc1,xc2,yc2,r1,r2){
  
  let distcentros = distPontoPonto(xc1, xc2, yc1, yc2);
  if (distcentros < Math.abs(r1-r2)){
    return "A circunferência é interna";
  }else if(distcentros == Math.abs(r1-r2)){
    return "A circunferência é tangente interna";
  }else if(distcentros < (r1+r2)){
    return "As circunferências são secantes";
  }else if(distcentros == (r1+r2)){
    return "A circunferência é tangente externa";
  }else if(distcentros >(r1+r2)){
    return "A circunferência é externa";
  }
}


// CINEMATICA (MU, MUV, MRU, MRUV e MCU) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------


function funcaoHorariaDeslocamento(S0, V, t){
  return (S0+V*t);

  /*
  S0 = localizacao inicial
  V = velocidade
  t = tempo
  */
}

function velMedia(Ds, Dt){
 return (Ds/Dt);

 /*
 Ds = variacao de espaco
 Dt = variacao de tempo
 */
}

function aceleracao(Dv, Dt){
  return (Dv/Dt);

  /*
  Dv = variacao da velocidade
  Dt = variacao de tempo
  */
}

function funcaoHorariaVel(V0, a, t){
  return (V0+a*t)

  /*
  V0 = velocidade inicial
  a = aceleracao
  t = tempo
  */
}

function funcaoHorariaPos(S0, V0, t, a){
  return (S0+V0*t+(1/2*a*Math.pow(t,2)));

  /*
  S0 = posicao inicial
  V0 = velocidade inicial
  t = tempo
  a = aceleracao
  */
}

function torricelli(V0, a, Ds){
  return Math.sqrt(Math.pow(V0,2)+(2*a*Ds));

  /*
  V0 = velocidade inicial
  a = aceleracao
  Ds = deslocamento
  */
}

function movimentoVertical(V0, g, h){
  return (Math.pow(V0,2)+(2*g*h));

  /*
  V0 = velocidade inicial
  g = aceleracao da gravidade
  h = altura
  */
}

function posicaoAngular(S,r){
  return (S/r);
  
  /*
  S = posicao
  r = raio da circunferencia
  */
}

function deslocAngular(Ds, r){
  return (Ds/r);

  /*
  Ds = variacao da posicao
  r = raio da circunferencia
  */
}

function velAngularMed(Da,Dt){
  return (Da/Dt);

  /*
  Da = deslocamento angular
  Dt = variacao do tempo
  */
}

function aceleracaoCentrip(V, r){
  return (Math.pow(V,2)/r);

  /*
  V = velocidade
  r = raio da circunferencia
  */
}

function aceleracaoAngMedia(w, Dt){
  return (w/Dt)

  /*
  w = velocidade angular
  Dt = variacao do tempo
  */
}

// DINAMICA (LEIS DE NEWTON, FORCA ATRITO E FORCA ELASTICA) ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function segundaLeiNewton(m,a){
  return (m*a);

  /*
  m = massa
  a = aceleracao
  */
}

function forcaPeso(m,g){
  return (m*g);

  /*
  m = massa
  g = aceeleracao da gravidadee
  */
}

function terceiraLeiNewton(F){
  return(F*-1);
  /*
  F = forca
  */
}

function forcaGravitacional(m1,m2,d){
  return(m1*m2*6.67*Math.pow(10,-11))/(Math.pow(d,2));

  /*
  m1 = massa do primeiro corpo
  m2 = massa do segundo corpo
  d = distancia entre os dois corpos
  */
}

function forcaAtrito(y,N){
  return(y*N);
  
  /*
  y = coeficiente de atrito
  N = forca normal
  */
}

function forcaElastica(K,x){
  return (K*x);

  /*
  K = constante elastica
  x = variacao de comprimento
  */
}

function energiaPotElastica(K,x){
  return (K*Math.pow(x,2))/2;

  /*
  K = constante elastica
  x = medida da deformacao do corpo elastico
  */
}


// OUTRAS FORMULAS DA MECANICA CLASSICA E FISICA (ENERGIA CINETICA, MOMENTO, TRABALHO, DENSIDADE E PRESSAO) -----------------------------------------------------------------------------------------------------------------------

function energiaCinetica(m,v){
  return (1/2*m*v*v);

  /*
  m = massa
  v = velocidade
  */
}

function momento(m,v){
  return (m*v);

  /*
  m = massa
  v = velocidade
  */
}

function trabalho(f,d){
  return (f*d);

  /*
  f = forca
  d = deslocamento
  */
}

function densidade(m,v){
  return (m/v);

  /*
  m = massa
  v = volume
  */
}

function pressao(f,a){
  return (f/a);

  /*
  f = forca
  a = area
  */
}

// TERMOLOGIA (CONVERSAO DE UNIDADES, DILATACAO TERMICA, QUANTIDADE DE CALOR E CAPACIDADE TERMICA) --------------------------------------------------------------------------------------------------------------------

function celsiusParaFahrenheit(c){
  return (c*(9/5))+32;
}

function fahrenheitParaCelsius(f){
  return (f-32)*(5/9);
}

function celsiusParaKelvin(c){
  return (c+273.15);
}

function kelvinParaCelsius(k){
  return (k-273.15);
}

function dilatacaoLinear(L0, a, Dt){
  return (L0*a*Dt);

  /*
  L0 = comprimento inicial
  a = coeficiente de dilatacao linear
  Dt = variacao de temperatura
  */
}

function dilatacaoSuperficial(A0, a, Dt){
  return (A0*2*a*Dt);

  /*
  A0 = area inicial
  a = coeficiente de dilatacao linear
  Dt = variacao de temperatura
  */
}

function dilatacaoVolumetrica(V0, a, Dt){
  return (V0*3*a*Dt);

  /*
  V0 = volume inicial
  a = coeficiente de dilatacao linear
  Dt = variacao de temperatura
  */
}

function qtdCalor(m, c, Dt){
  return (m*c*Dt);
  /*
  m = massa
  c = calor especifico
  Dt = variacao de temperatura
  */
}

function capacidadeTermica(m,c){
  return (m*c);

  /*
  m = massa
  c = calor especifico
  */
}

// OPTICA (POTENCIA FOCAL, EQUACAO DE GAUSS, ASSOCIACAO DE ESPELHOS PLANOS, AMPLIACAO E INDICE DE REFRACAO) -----------------------------------------------------------------------------------------------------------

function potenciaFocal(f){
  return (1/f);
  // f = distancia focal
}

function equacaoGauss(Di, Do){
  return 1/(Di+Do);

  /*
  Di = distancia da imagem
  Do = distancia do objeto
  */
}

function assocEspelhosPlanos(a){
  return (360/a)-1;

  // a = angulo de abertura entre os espelhos
  // resultado sera o numero de imagens entre esses 2 espelhos
}

function ampliacao(f, Do){
  return(f/(f-Do));

  /*
  f = distancia focal
  Do = distancia do objeto
  */
}

function indiceRefracao(v){
  return (299792458/v);

  /*
  v = velocidade da luz no meio em m/s
  */
}

// ONDULATORIA (PERIODO E VELOCIDADE DE PROPAGACAO)

function velPropOndaLinear(T, m, l){
  return Math.sqrt(T/(m/l));

  /*
  T = tracao
  m = massa
  l = comprimento
  */
}

function periodo(f){
  return 1/f;

  // f = frequencia
}

function velPropagacao(l,f){
  return l*f;

  /*
  l = comprimento de onda
  f = frequencia
  */
}

// ELETROSTATICA (QUANTIDADE DE CARGA ELETRICA, LEI DE COULOMB E CAMPO ELETRICO) ---------------------------------------------------------------------------------------------------------------------------------------

function qtdCargaEletrica(n){
  return n*(1.6*Math.pow(10, -19));

  // n = numero de eletrons em falta ou sobrando  
}

function leiCoulomb(q1, q2, d){
  return (q1*q2*(9*Math.pow(10,9)))/(Math.pow(d,2));

  /*
  q1 = carga 1
  q2 = carga 2
  d = distancia entre as cargas
  */
}

function campoEletrico(F, q){

  return (F/q);

  /*
  F = forca
  q = carga
  */
}

// ELETRODINAMICA (CORRENTE ELETRICA, PRIMEIRA E SEGUNDA LEIS DE OHM, POTENCIA ELETRICA, ENERGIA ELETRICA, RESISTENCIA EQUIVALENTE EM SERIE E PARALELO)

function correnteEletrica(Q, t){
  return (Q/t);

  /*
  Q = quantidade de carga eletrica
  t = intervalo de tempo

  */
}

function primeiraLeiOhm(R, i){
  return (R*i);

  /*
  R = resistencia
  i = corrente
  */
}

function segundaLeiOhm(p,l,A){
return ((p*l)/A);

/*
p = resistividade
l = comprimento do fio
A = area transversal do fio (bitola)
*/
}

function potenciaEletrica(i, U){
  return (i*U);
  /*
  i = corrente
  U = tensao
  */
}

function potenciaEletricaSemTensao(R, i){
  return (R*i*i);

  /*
  R = resistencia
  i = corrente
  */
}

function potenciaEletricaSemCorrente(U,R){

  return ((U*U)/R);

  /*
  U = tensao
  R = resistencia
  */
}

function energiaEletrica(P,t){
  return P*t;

  /*
  P = potencia
  t = intervalo de tempo
  */
}

function ReqSerie(){
  

  var Req = parseInt(0);
  var limite = prompt("Digite a quantidade de resistores");
  
  for(var i = 0; i < limite; i++){
     var R = parseInt((prompt("Digite o valor do resistor")));
   
     var Req = Req + R;

  }

  return ("O valor da resistência equivalente é de "+Req+" Ohms");
}


function ReqParalelo(){
  var Req = parseInt(0);
  var limite = prompt("Digite a quantidade de resistores");
  
  for(var i = 0; i < limite; i++){
     var R = parseInt((prompt("Digite o valor do resistor")));
   
     var Req = Req + (1/R);

  }

  return ("O valor da resistência equivalente é de "+Math.pow(Req,-1))+" Ohms";

}

// EDUCACAO FISICA (IMC, TORNEIO MATA-MATA E TORNEIO ROUND-ROBIN)

function imc(p,a){

  return (p/(a*a));

  /*
  p = peso
  a = altura
  */
}

function mataMata(n){

  return (n-1);

  /*
  n = numero de times participantes
  */
}

function roundRobin(n){

  return combinacao(n,2);

  /*
  n = numero de times participantes
  */
}

// FAVOR SUBSTITUIR O CODIGO 'return' PELO CODIGO 'res.innerHTML =' QUANDO FOR NECESSARIO TESTAR
// ALEM DISSO, ACRESCENTAR AS FUNCOES AS DEFINICOES DAS VARIAVEIS, EXEMPLO: 'var n = parseFloat(document.getElementById("n").value);'
