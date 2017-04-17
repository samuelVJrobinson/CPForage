#Function to find  S (scaling term), given Load, Forager number, Max intake rate, travel time, and nectar production
scalFun=function(mu=NULL,l=NULL,NumFls=NULL,L_i=NULL,L_max_i=NULL,n_i=NULL,h_i=NULL,p_i=NULL,f_i=NULL,d_i=NULL,v_i=NULL,beta_i=NULL,H_i=NULL,patchLev=F){
  #Argument checking
  singleArgs=list(mu=mu,l=l,NumFls=NumFls) #Patch arguments (mu,l,NumFls)
  singArgLen=sapply(singleArgs,length) #Length of patch arguments
  if(sum(singArgLen>1)>0){ #If patch arguments are longer than 1
    stop(paste(names(singArgLen)[singArgLen>1],'longer than 1. Patch arguments must be a scalar.'))
  }
  if(singleArgs$mu==0|singleArgs$l==0|singleArgs$NumFls==0) return(NA) #If production is 0, competition isn't defined.
  #Forager (nest-level) arguments (Load,MaxLoad,number,Handling Time,Licking Rate,b/w flower flight time,distance,beta,hive time)
  Args=list(L_i=L_i,L_max_i=L_max_i,n_i=n_i,h_i=h_i,p_i=p_i,f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,H_i=H_i)
  argLen=sapply(Args,length) #Length of forager arguments
  argSign=sapply(c(singleArgs,Args),function(x) sum(x<0))>0 #Are any arguments less than 0?
  if(sum(argSign)>0){ #If any arguments are <0
    stop(c(singleArgs,Args)[argSign],'< 0')
  } else if((sum(argLen==0)+sum(singArgLen==0))>0) {#If there are any zero-length arguments (i.e. not provided)
    stop(paste(names(singArgLen)[singArgLen==0],names(argLen)[argLen==0],'arguments not provided'))
  } else if(length(unique(argLen))>2) { #If there are >2 lengths of argument
    stop(paste(length(unique(argLen)),'different vector lengths provided to scaleFun. Requires either 2 or 1.'))
  } else if(length(unique(argLen))==2){ #If there are 2 lengths of argument
    if(sum(unique(argLen)==1)==0){ #If neither argument is 1 long
      stop('2 different vector lengths provided, but neither is 1 long.')
    }
    lmax=max(unique(argLen))
    for(i in names(argLen)[argLen==1]) assign(i,rep(get(i),lmax))
  }
  #Function for calculating patch usage or flower depletion
  usage=function(S,mu,l,NumFls,L_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,L_max_i,patchLev){
    if(patchLev){ #Patch-level version (uL/s)
    return(sum((L_i*n_i)/((L_i*(h_i+S*l*p_i+f_i)/S*l)+(d_i/v_i)*((2-beta_i*L_i/L_max_i)/(1-beta_i*L_i/L_max_i))+H_i)))
    } else { #Flower-level version (uL per flower)
    return(sum(l/((L_i*n_i)/(mu*NumFls*S*((d_i*(2*L_max_i-L_i*beta_i))/(v_i*(L_max_i-L_i*beta_i))+(L_i*(S*l*p_i+h_i+f_i))/(S*l)+H_i))+1)))
    }
  }
  #Function for finding S by optimization
  usageS=function(S,mu,l,NumFls,L_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,L_max_i,patchLev){
    if(patchLev){ #Patch-level version
      return(abs(usage(S,L_i,n_i,h_i,l,p_i,f_i,d_i,v_i,beta_i,L_max_i)-(mu*NumFls))) #Usage must = mu*NumFls
    } else { #Flower-level version
      return(abs(usage(S,mu,l,NumFls,L_i,n_i,h_i,p_i,f_i,d_i,v_i,beta_i,L_max_i)-S*l)) #Mean nectar must = S*l
    }
  }
  lmax=max(unique(argLen)) #Length of arguments (number of cells to process)
  if(patchLev){ #Patch-level S-value
    maxUse=usage(1,L_i,n_i,h_i,l,p_i,f_i,d_i,v_i,beta_i,L_max_i,patchLev) #Patch-level usage when S=1
    if((mu*NumFls)>=maxUse){ #If patch production >= max usage
      return(1) #No competition
    } else { #If max usage > patch production
      if(lmax==1){ #If only 1 nest is present
        #Original version
        # S=(L_i*(h_i+f_i))/(l*((L_i*n_i/mu)-(d_i*(2-beta_i*L_i/L_max_i)/v_i*(1-beta_i*L_i/L_max_i))-H_i-L_i*p_i))
        #Version from Maxima
        S=-(L_i*(L_max_i-beta_i*L_i)*(h_i+f_i)*(mu*NumFls)*v_i)/(l*(L_i*L_max_i*(mu*NumFls)*p_i*v_i-beta_i*L_i^2*(mu*NumFls)*p_i*v_i-L_i*L_max_i*n_i*v_i+beta_i*L_i^2*n_i*v_i+H_i*L_max_i*(mu*NumFls)*v_i-beta_i*H_i*L_i*(mu*NumFls)*v_i+2*L_max_i*d_i*(mu*NumFls)-beta_i*L_i*d_i*(mu*NumFls)))
      } else { #If there are multiple nests
        S=optimize(usageS,interval=c(0,1),L_i=L_i,n_i=n_i,h_i=h_i,l=l,p_i=p_i,f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,L_max_i=L_max_i,mu=mu,NumFls=NumFls)$min #Finds S by optimization
      }
    }
  } else { #Flower-level S-value using Possingham 1988
    #In this version, flower-level depletion begins almost immediately
    if(lmax==1){ #If only 1 nest is present
      #Version from Maxima -
      S=-(sqrt(((mu^2*L_i^4*NumFls^2*l^2*p_i^2+(-2*mu*L_i^4*NumFls*l^2*n_i+2*mu^2*H_i*L_i^3*NumFls^2*l^2+
        (2*mu^2*L_i^4*NumFls^2*h_i+2*mu^2*L_i^4*NumFls^2*f_i)*l)*p_i+L_i^4*l^2*n_i^2+((2*mu*L_i^4*NumFls*h_i+
        2*mu*L_i^4*NumFls*f_i)*l-2*mu*H_i*L_i^3*NumFls*l^2)*n_i+mu^2*H_i^2*L_i^2*NumFls^2*l^2+(2*mu^2*H_i*L_i^3*
        NumFls^2*h_i+2*mu^2*H_i*L_i^3*NumFls^2*f_i)*l+mu^2*L_i^4*NumFls^2*h_i^2+2*mu^2*L_i^4*NumFls^2*f_i*h_i+mu^2*
        L_i^4*NumFls^2*f_i^2)*v_i^2+(2*mu^2*L_i^3*NumFls^2*d_i*l^2*p_i-2*mu*L_i^3*NumFls*d_i*l^2*n_i+2*mu^2*H_i*
        L_i^2*NumFls^2*d_i*l^2+(2*mu^2*L_i^3*NumFls^2*d_i*h_i+2*mu^2*L_i^3*NumFls^2*d_i*f_i)*l)*v_i+mu^2*L_i^2*
        NumFls^2*d_i^2*l^2)*beta_i^2+((-2*mu^2*L_i^3*L_max_i*NumFls^2*l^2*p_i^2+(4*mu*L_i^3*L_max_i*NumFls*l^2*n_i-4*
        mu^2*H_i*L_i^2*L_max_i*NumFls^2*l^2+(-4*mu^2*L_i^3*L_max_i*NumFls^2*h_i-4*mu^2*L_i^3*L_max_i*NumFls^2*f_i)*l)*
        p_i-2*L_i^3*L_max_i*l^2*n_i^2+(4*mu*H_i*L_i^2*L_max_i*NumFls*l^2+(-4*mu*L_i^3*L_max_i*NumFls*h_i-4*mu*L_i^3*
        L_max_i*NumFls*f_i)*l)*n_i-2*mu^2*H_i^2*L_i*L_max_i*NumFls^2*l^2+(-4*mu^2*H_i*L_i^2*L_max_i*NumFls^2*h_i-4*mu^2*
        H_i*L_i^2*L_max_i*NumFls^2*f_i)*l-2*mu^2*L_i^3*L_max_i*NumFls^2*h_i^2-4*mu^2*L_i^3*L_max_i*NumFls^2*f_i*h_i-2*
        mu^2*L_i^3*L_max_i*NumFls^2*f_i^2)*v_i^2+(-6*mu^2*L_i^2*L_max_i*NumFls^2*d_i*l^2*p_i+6*mu*L_i^2*L_max_i*NumFls*
        d_i*l^2*n_i-6*mu^2*H_i*L_i*L_max_i*NumFls^2*d_i*l^2+(-6*mu^2*L_i^2*L_max_i*NumFls^2*d_i*h_i-6*mu^2*L_i^2*
        L_max_i*NumFls^2*d_i*f_i)*l)*v_i-4*mu^2*L_i*L_max_i*NumFls^2*d_i^2*l^2)*beta_i+(mu^2*L_i^2*L_max_i^2*NumFls^2*
        l^2*p_i^2+(-2*mu*L_i^2*L_max_i^2*NumFls*l^2*n_i+2*mu^2*H_i*L_i*L_max_i^2*NumFls^2*l^2+(2*mu^2*L_i^2*
        L_max_i^2*NumFls^2*h_i+2*mu^2*L_i^2*L_max_i^2*NumFls^2*f_i)*l)*p_i+L_i^2*L_max_i^2*l^2*n_i^2+((2*mu*L_i^2*
        L_max_i^2*NumFls*h_i+2*mu*L_i^2*L_max_i^2*NumFls*f_i)*l-2*mu*H_i*L_i*L_max_i^2*NumFls*l^2)*n_i+mu^2*H_i^2*
        L_max_i^2*NumFls^2*l^2+(2*mu^2*H_i*L_i*L_max_i^2*NumFls^2*h_i+2*mu^2*H_i*L_i*L_max_i^2*NumFls^2*f_i)*l+mu^2*
        L_i^2*L_max_i^2*NumFls^2*h_i^2+2*mu^2*L_i^2*L_max_i^2*NumFls^2*f_i*h_i+mu^2*L_i^2*L_max_i^2*NumFls^2*f_i^2)*
        v_i^2+(4*mu^2*L_i*L_max_i^2*NumFls^2*d_i*l^2*p_i-4*mu*L_i*L_max_i^2*NumFls*d_i*l^2*n_i+4*mu^2*H_i*L_max_i^2*
        NumFls^2*d_i*l^2+(4*mu^2*L_i*L_max_i^2*NumFls^2*d_i*h_i+4*mu^2*L_i*L_max_i^2*NumFls^2*d_i*f_i)*l)*v_i+4*mu^2*
        L_max_i^2*NumFls^2*d_i^2*l^2)+((-mu*L_i^2*NumFls*l*p_i+L_i^2*l*n_i-mu*H_i*L_i*NumFls*l+mu*L_i^2*NumFls*h_i+mu*
        L_i^2*NumFls*f_i)*v_i-mu*L_i*NumFls*d_i*l)*beta_i+(mu*L_i*L_max_i*NumFls*l*p_i-L_i*L_max_i*l*n_i+mu*H_i*
        L_max_i*NumFls*l-mu*L_i*L_max_i*NumFls*h_i-mu*L_i*L_max_i*NumFls*f_i)*v_i+2*mu*L_max_i*NumFls*d_i*l)/(((2*mu*
        L_i^2*NumFls*l*p_i+2*mu*H_i*L_i*NumFls*l)*v_i+2*mu*L_i*NumFls*d_i*l)*beta_i+(-2*mu*L_i*L_max_i*NumFls*l*p_i-2*
        mu*H_i*L_max_i*NumFls*l)*v_i-4*mu*L_max_i*NumFls*d_i*l)
    } else { #If there are multiple nests
      S=optimize(usageS,interval=c(0,1),L_i=L_i,n_i=n_i,h_i=h_i,l=l,p_i=p_i,f_i=f_i,d_i=d_i,v_i=v_i,beta_i=beta_i,L_max_i=L_max_i,mu=mu)$min #Find S by optimization
    }
  }
  return(S)
}
