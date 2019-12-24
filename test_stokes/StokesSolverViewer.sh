VAR=$1
NUM=$2

echo $0".Variable" : $VAR
echo $0".TimeStep" : $NUM

BASENAMESOL=integrator.$VAR.`printf "%.5d\n" $NUM`

if [ $VAR == "p" ];then
    echo $0".DisplayPressure"
    ln -fs integrator.P1.00000.mesh $BASENAMESOL.mesh
fi

if [ $VAR == "u" ];then
    echo $0".DisplayVelocity"
    ln -fs integrator.P2.00000.mesh $BASENAMESOL.mesh
fi

MnsViewer $BASENAMESOL
\rm $BASENAMESOL.mesh


