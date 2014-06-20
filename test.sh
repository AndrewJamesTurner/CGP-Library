#!/bin/sh

make gettingStarted
make createDataSet
make manipulatingChromosomes
make customNodeFunction
make customFitnessFunction
make customSelectionScheme
make customReproductionScheme
make averageBehaviour
make neuroEvolution
make printChromoEqu
make customES
make visualization

echo "\n\n ** gettingStarted ** \n\n"
valgrind --leak-check=yes --error-exitcode=1 ./gettingStarted > /dev/null

if [ $? -eq 1 ]; then
	echo gettingStarted: memory leak...
	exit
else
	echo gettingStarted: PASS
fi



echo "\n\n ** createDataSet ** \n\n"
valgrind --leak-check=yes --error-exitcode=1 ./createDataSet > /dev/null

if [ $? -eq 1 ]; then
	echo createDataSet: memory leak...
	exit
else
	echo createDataSet: PASS
fi


echo "\n\n ** manipulatingChromosomes ** \n\n"
valgrind --leak-check=yes --error-exitcode=1 ./manipulatingChromosomes > /dev/null

if [ $? -eq 1 ]; then
	echo manipulatingChromosomes: memory leak...
	exit
else
	echo manipulatingChromosomes: PASS
fi


echo "\n\n ** customNodeFunction ** \n\n"
valgrind --leak-check=yes --error-exitcode=1 ./customNodeFunction > /dev/null

if [ $? -eq 1 ]; then
	echo customNodeFunction: memory leak...
	exit
else
	echo customNodeFunction: PASS
fi


echo "\n\n ** customFitnessFunction ** \n\n"
valgrind --leak-check=yes --error-exitcode=1 ./customFitnessFunction > /dev/null

if [ $? -eq 1 ]; then
	echo customFitnessFunction: memory leak...
	exit
else
	echo customFitnessFunction: PASS
fi

echo "\n\n ** customSelectionScheme ** \n\n"
valgrind --leak-check=yes --error-exitcode=1 --error-exitcode=1 ./customSelectionScheme > /dev/null

if [ $? -eq 1 ]; then
	echo customSelectionScheme: memory leak...
	exit
else
	echo customSelectionScheme: PASS
fi

echo "\n\n ** customReproductionScheme ** \n\n"
valgrind --leak-check=yes --error-exitcode=1 ./customReproductionScheme > /dev/null

if [ $? -eq 1 ]; then
	echo customReproductionScheme: memory leak...
	exit
else
	echo customReproductionScheme: PASS
fi

echo "\n\n ** averageBehaviour ** \n\n"
valgrind --leak-check=yes --error-exitcode=1 ./averageBehaviour > /dev/null

if [ $? -eq 1 ]; then
	echo averageBehaviour: memory leak...
	exit
else
	echo averageBehaviour: PASS
fi

echo "\n\n ** neuroEvolution ** \n\n"
valgrind --leak-check=yes --error-exitcode=1 ./neuroEvolution > /dev/null

if [ $? -eq 1 ]; then
	echo neuroEvolution: memory leak...
	exit
else
	echo neuroEvolution: PASS
fi

echo "\n\n ** printChromoEqu ** \n\n"
valgrind --leak-check=yes --error-exitcode=1 ./printChromoEqu > /dev/null

if [ $? -eq 1 ]; then
	echo printChromoEqu: memory leak...
	exit
else
	echo printChromoEqu: PASS
fi

rm tmp.tex

echo "\n\n ** customES ** \n\n"
valgrind --leak-check=yes --error-exitcode=1 ./customES > /dev/null

if [ $? -eq 1 ]; then
	echo customES: memory leak...
	exit
else
	echo customES: PASS
fi


echo "\n\n ** visualization ** \n\n"
valgrind --leak-check=yes --error-exitcode=1 ./visualization > /dev/null

if [ $? -eq 1 ]; then
	echo visualization: memory leak...
	exit
else
	echo visualization: PASS
fi

echo "\n\n Test Complete :) \n\n"


rm chromo.*

make clean
