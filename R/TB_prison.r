TB_prison <- function (t, state, parameters) {
  with(as.list(c(state, parameters)),             
       {
        N       <- U + L + I + R                 # Total population
        Np      <- Up + Lp + Ip + Rp            # Total population
        Nx      <- Ux + Lx + Ix + Rx            # Total population
        births <- (I+Ip+Ix) * mu_tb + (N+Np+Nx) * mu  # Births (for stable population)
        lambda   <- beta * (I+Ix)/(N+Nx)        # Force of Infection coomunirt
        lambda_p <- beta_p * Ip/Np              # Force of Infection prisons

       
   
         ##### Community transmission
         # Uninfected 
         dU <- births - U * (lambda + mu + r_incar)                              
         
         # Latent 
         dL <- U * lambda * (1-fast) + R * (lambda * (1-fast) * imm) - L * (mu + break_in + r_incar) 
         
         # Active TB
         dI <- U * lambda * fast + R * (lambda * fast * imm) +  L * break_in + R * relapse -
               I * (mu + mu_tb + selfcure + r_incar)
         
         # Recovered
         dR <- I * selfcure - R * (lambda * imm + relapse + mu + r_incar)   
         
         
         ###### Ex-prison transmission
         # Uninfected 
         dUx <- Up * r_release -  Ux * (lambda + mu + r_reincar)                              
         
         # Latent 
         dLx <- Ux * lambda * (1-fast) + Rx * (lambda * (1-fast) * imm) + 
           Lp * r_release - Lx * (mu + break_in_x + r_reincar) 
         
         # Active TB
         dIx <- Ux * lambda * fast + Rx * (lambda * fast * imm) + 
           Lx * break_in_x + Rx * relapse + Ip * r_release  - 
           Ix * (mu + mu_tb + selfcure_x + r_reincar)
         
         # Recovered
         dRx <- Ix * selfcure_x + Rp * r_release - 
           Rx * (lambda * imm + relapse + mu + r_reincar)      
         
 
         ## Prison Transmission

          # Uninfected 
         dUp <- U * r_incar + Ux * r_reincar -  Up * (lambda_p + mu + r_release)                              
         
         # Latent 
         dLp <- Up * lambda_p * (1-fast_p) + Rp * (lambda_p * (1-fast_p) * imm) + 
           L * r_incar + Lx * r_reincar - Lp * (mu + break_in_p + r_release) 
         
         # Active TB
         dIp <- Up * lambda_p * fast_p + Rp * (lambda_p * fast_p * imm) +  
           Lp * break_in_p + Rp * relapse + 
           I * r_incar + Ix * r_reincar - Ip * (mu + mu_tb + selfcure_p + r_release)
         
         # Recovered
         dRp <- Ip * selfcure_p + R * r_incar + Rx * r_reincar - 
           Rp * (lambda_p * imm + relapse + mu + r_release)      

         
         # Model outcomes

         # Community 
         dIncidence <- U * lambda * fast + R * (lambda * fast * imm) +  L * break_in + R * relapse
                  
         dIrecent   <- U * (lambda * fast) + R * (lambda * fast * imm)

         dIremote   <-L * break_in + R * relapse
         
         # Ex Prison 
         dIncidence_x <- Ux * (lambda * fast) + Rx * (lambda * fast * imm) + Lx * break_in_x + Rx * relapse 
         
         dIrecent_x   <-  Ux * lambda * fast + Rx * (lambda * fast * imm) + Lx * break_in_x + Rx * relapse + Ip * r_release
         
         dIremote_x   <- Lx * break_in_x + Rx * relapse
         

        # Prison 
         dIncidence_p <-  Up * lambda_p * fast_p + Rp * (lambda_p * fast_p * imm) +  
           Lp * break_in_p + Rp * relapse + I * r_incar + Ix * r_reincar 
                  
         dIrecent_p   <- Up * (lambda_p * fast_p) + Rp * (lambda_p * fast_p * imm)

         dIremote_p   <-Lp * break_in_p + Rp * relapse
         

        
         # ::::::::::::::::::::::::::::::::::::::::::::
           
         # wrap-up 
        dx <- c(dU, dL, dI, dR,
                dUx, dLx, dIx, dRx,
                dUp, dLp, dIp, dRp,
                dIncidence, dIrecent, dIremote,
                dIncidence_x, dIrecent_x, dIremote_x,
                dIncidence_p, dIrecent_p, dIremote_p)
         list(dx)
       }
  )
}
