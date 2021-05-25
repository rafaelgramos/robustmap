# Robustmap

For a version in English, check [README.md](https://github.com/rafaelgramos/robustmap/blob/master/README.md).

## Visão Geral

Este pacote para a linguagem R implementa rotinas para mapear conjuntos de pontos geocodificados (e.g. crimes, casos de uma doença) via contagem de quadrantes e densidade *kernel* levando-se em consideração a robustez e a uniformidade espacial interna das contagens agregadas. Em outras palavras, quadrantes muito grandes (ou uma largura de banda muito longa) podem levar a uma falta de detalhes e à Falácia Ecológica; por outro lado, quadrantes muito pequenos (ou uma largura de banda muito curta) podem gerar contagens e densidades pouco confiáveis (muito vulneráveis a ruído/variações aleatórias).

A rotina **robust.quadcount** estima um tamanho de quadrante que equilibre robutez e uniformidade interna (dado um conjunto de pontos georrefenciados), retornando uma mapa de contagem por quadrantes 'otimizado', além de outras informações como: quais a granularidades (tamanho de quadrante) que foram testadas, a robustez e uniformidade interna para cada granularidade, etc. Esta metodologia é baseada no estudo de [Ramos, et al. (2020)](https://link.springer.com/article/10.1007/s10940-020-09474-6)

A rotina **robust.density** está ainda sendo desenvolvida e deverá estimar a densidade de um padrão de pontos georreferenciados usando um largura de banda variável, em que a largura de banda é a maxima possível tal que o subconjunto de pontos contidados nela possua Aleatoriedade Espacial Total. Desta forma, para qualquer pixel, o número de pontos considerados na estimação da densidade é maximado (melhorand a robustez) ao mesmo tempo evitando a Falácia Ecológica de tomar a densidade média numa área em que a densidade varia.

## Instalação

No Mac (ainda não testado para Linux ou Windows, mas deve funcionar também), primeiro instale o pacate **devtools** rodando a seguinte linha em um console do R :

	install.packages("devtools")

Então, ainda no console do R, instale of pacote **robustmap** a partir do GitHub rodando a liha:

	devtools::install_github("rafaelgramos/robustmap")

Após isto, você deverá ser capaz de carregar o pacote normalmente.