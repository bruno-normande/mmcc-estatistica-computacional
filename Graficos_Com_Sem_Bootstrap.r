plot(base_U,Distribuicao_U,type='l',col='blue',xlab="Distribuição U",ylab="Valores")
lines(base_U,Distribuicao_U_B,type='l',col='red')
legend(5,0.05,c("Com Bootstrap","Sem Bootstrap"),lty=c(1,1,1), lwd=c(0.5,0.5,0.5),col=c('red','blue'))

plot(base_BN,Distribuicao_BN,type='l',col='blue',xlab="Distribuição BN",ylab="Valores")
lines(base_BN,Distribuicao_BN_B,type='l',col='red')
legend(1,0.12,c("Com Bootstrap","Sem Bootstrap"),lty=c(1,1,1), lwd=c(0.5,0.5,0.5),col=c('red','blue'))

plot(base_Exp,Distribuicao_Exp,type='l',col='blue',xlab="Distribuição Exp",ylab="Valores")
lines(base_Exp,Distribuicao_Exp_B,type='l',col='red')
legend(1,-0.2,c("Com Bootstrap","Sem Bootstrap"),lty=c(1,1,1), lwd=c(0.5,0.5,0.5),col=c('red','blue'))
