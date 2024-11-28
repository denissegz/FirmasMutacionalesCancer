# Análisis de Firmas Mutacionales en Diferentes Tipos de Cáncer

El análisis realizado en este estudio tuvo como objetivo principal identificar las **firmas mutacionales características** de diferentes tipos de cáncer para mejorar el diagnóstico y la clasificación de estos tumores. Para ello, se emplearon diversos clasificadores, en particular el modelo **Random Forest**, utilizando mutaciones que se repiten en al menos cinco pacientes para cada tipo de cáncer. Este enfoque permitió evaluar la eficacia de los modelos para diferenciar entre el cáncer de interés (clase 1) y otros tipos de cáncer (clase 0). A continuación, se presentan las conclusiones en función de los resultados obtenidos:

## 1. Identificación de Firmas Mutacionales para Diferentes Tipos de Cáncer

Se analizaron 13 tipos de cáncer, y se encontró que las **firmas mutacionales** de ciertos cánceres, como los de **tiroides**, **hígado** y **riñón**, fueron altamente distinguibles mediante el modelo Random Forest. En estos casos, el número de mutaciones específicas y la cantidad de pacientes para cada tipo de cáncer fueron suficientes para entrenar modelos robustos que alcanzaron buenas métricas en términos de **precisión**, **sensibilidad** y **especificidad**.

- El **cáncer de tiroides** mostró una gran capacidad para ser identificado con una alta **sensibilidad (99.7%)** y **especificidad (85.4%)**, con un **área bajo la curva ROC (AUC)** de **0.986**.
- El **cáncer de hígado**, a pesar de tener una baja cantidad de mutaciones (**8 mutaciones en 53 pacientes**), demostró también un alto rendimiento, con una **AUC de 0.937** y una **sensibilidad del 99.6%**.
- El **cáncer de riñón**, con solo 3 mutaciones en 14 pacientes, presentó un modelo perfecto (**AUC de 1.0**) en el conjunto de entrenamiento y prueba.

## 2. Desafíos con Cánceres de Baja Frecuencia de Mutaciones

Los cánceres con menor número de mutaciones y/o pacientes, como el **cáncer de glándula adrenal**, presentaron resultados menos fiables, con **especificidad baja** en el modelo, lo que dificultó la separación de las mutaciones específicas de estos cánceres de las del resto de los cánceres. 

Además, los modelos mostraron dificultades de generalización cuando los cánceres no presentaron mutaciones claramente distintivas o cuando había una alta superposición en el espacio de mutaciones. Esto se observó particularmente en cánceres como el **cáncer de esófago** y el **cáncer de estómago**, que, a pesar de contar con una mayor cantidad de pacientes, no mostraron un rendimiento óptimo en los clasificadores.

## 3. Validación de los Modelos

El conjunto de validación, que incluyó un total de **8133 pacientes** y **4569 mutaciones**, mostró que los modelos entrenados para **cáncer de tiroides**, **cáncer de hígado** y **cáncer de riñón** continuaron mostrando buenos resultados, aunque con variaciones en las métricas:

- Para el **cáncer de tiroides**, el modelo en validación mostró un desempeño aceptable con una **precisión de 97.6%**, **sensibilidad de 98.96%** y **especificidad de 70.2%**, con **AUC de 0.89**. Ajustando el umbral de decisión, se mejoró la especificidad a **78.3%**, sin sacrificar demasiado la sensibilidad.
- Sin embargo, el modelo para el **cáncer de hígado** mostró un descenso en la **especificidad (21.2%)**, lo que sugiere que el modelo tenía dificultad para diferenciar entre las mutaciones del hígado y las de otros cánceres, con un **AUC de 0.60**, lo que indica un rendimiento moderado.
- Para el **cáncer de riñón**, aunque la **sensibilidad alcanzó 100%**, la **especificidad fue nula**, lo que sugiere que el modelo pudo predecir correctamente los casos positivos pero falló al identificar los negativos, con un **AUC de 0.5**, similar a la aleatoriedad.

## 4. Importancia de las Mutaciones en la Clasificación

En todos los modelos, las **mutaciones específicas** para cada cáncer mostraron ser las variables más importantes para la clasificación. Las **4 mutaciones** para el cáncer de tiroides, las **8** para el cáncer de hígado y las **3** para el cáncer de riñón fueron identificadas como las más relevantes. Esto respalda la idea de que algunas mutaciones son fuertemente indicativas de la presencia de ciertos tipos de cáncer y pueden servir como **marcadores diagnósticos confiables**.

## Conclusiones Finales

Este estudio demuestra que es posible utilizar **firmas mutacionales específicas** para distinguir diferentes tipos de cáncer, con un enfoque efectivo mediante el uso de modelos de **aprendizaje automático**. Los resultados obtenidos para cánceres como el tiroides, hígado y riñón son prometedores, destacando la utilidad de identificar mutaciones características en estos cánceres para su diagnóstico temprano.

Sin embargo, también se identificaron **desafíos importantes**, como la baja especificidad en ciertos cánceres y la dificultad para manejar casos con un número limitado de mutaciones y pacientes. La mejora de estos modelos y la inclusión de más datos y mutaciones podrían ayudar a superar estos retos, optimizando la capacidad de los clasificadores para predecir con mayor precisión la presencia de cáncer en nuevos pacientes.
