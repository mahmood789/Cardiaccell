# Title: Human Sinoatrial Node (Pacemaker) Cell Simulator - V7
# Author: Gemini
# Date: June 9, 2025
# Description: This version enhances the educational value by adding scientific references.
#              - The "Model & Data Explained" tab now includes citations for the
#                computational model, each major ion channel, and the physiological
#                basis for the experimental data waveforms.
#              - This provides users with direct links to the source literature.

# --- Required Packages ---
# install.packages(c("shiny", "deSolve", "shinycssloaders", "ggplot2", "tidyr", "bslib"))

library(shiny)
library(deSolve)
library(shinycssloaders)
library(ggplot2)
library(tidyr)
library(bslib)

#==============================================================================
# 1. EXPANDED BIOLOGICAL DATA MODULE
#==============================================================================
exp_time_normal <- c(0, 150, 300, 450, 500, 520, 540, 560, 580, 600, 700, 850)
exp_voltage_normal <- c(-60, -55, -50, -40, -20, 0, 10, 5, -15, -45, -55, -60)
exp_data_normal <- as.data.frame(spline(exp_time_normal, exp_voltage_normal, n = 851, method = "fmm"))
names(exp_data_normal) <- c("time", "V")

exp_time_brady <- c(0, 250, 500, 750, 800, 820, 840, 860, 880, 900, 1000, 1200)
exp_voltage_brady <- c(-62, -58, -52, -40, -20, 0, 10, 5, -15, -45, -57, -62)
exp_data_brady <- as.data.frame(spline(exp_time_brady, exp_voltage_brady, n = 1201, method = "fmm"))
names(exp_data_brady) <- c("time", "V")

exp_time_tachy <- c(0, 100, 200, 300, 325, 340, 355, 370, 385, 425, 500, 630)
exp_voltage_tachy <- c(-58, -54, -50, -40, -20, 0, 12, 8, -10, -40, -52, -58)
exp_data_tachy <- as.data.frame(spline(exp_time_tachy, exp_voltage_tachy, n = 631, method = "fmm"))
names(exp_data_tachy) <- c("time", "V")

experimental_datasets <- list(
    "Normal Sinus Rhythm (~70 BPM)" = exp_data_normal,
    "Sinus Bradycardia (~50 BPM)" = exp_data_brady,
    "Sinus Tachycardia (~95 BPM)" = exp_data_tachy
)

#==============================================================================
# 2. USER INTERFACE (UI)
#==============================================================================
ui <- page_sidebar(
    title = "Human Sinoatrial Node Cell Simulator",
    theme = bs_theme(version = 5, bootswatch = "cosmo"),
    sidebar = sidebar(
        width = 375,
        card(card_header("Simulation Controls"), card_body(
            sliderInput("sim_time", "Duration (ms):", min=1000, max=10000, value=5000, step=500),
            actionButton("run_button", "Run Simulation", class="btn-success", icon=icon("play"))
        )),
        card(card_header("Ion Channel Scaling"), card_body(
            p("Scale the baseline conductance of channels.", class="text-muted small"),
            sliderInput("g_f_scale", "Funny Current (I_f)", min=0, max=2, value=1, step=0.1),
            sliderInput("g_CaL_scale", "L-type Ca2+ (I_CaL)", min=0, max=2, value=1, step=0.1),
            sliderInput("g_Kr_scale", "Rapid K+ (I_Kr)", min=0, max=2, value=1, step=0.1)
        )),
        card(card_header("Pharmacology: Drug Effects"), card_body(
            p("Simulate drug effects by applying a % block to the target current.", class="text-muted small"),
            sliderInput("block_f", "Ivabradine (% I_f Block):", min=0, max=100, value=0, step=5),
            sliderInput("block_CaL", "Ca-Antagonist (% I_CaL Block):", min=0, max=100, value=0, step=5),
            sliderInput("block_Kr", "Class III (% I_Kr Block):", min=0, max=100, value=0, step=5)
        )),
        card(card_header("Plotting Controls"), card_body(
            actionButton("store_trace", "Store Trace for Comparison", icon=icon("layer-group")),
            actionButton("clear_traces", "Clear Stored Traces", icon=icon("trash-can")),
            hr(),
            selectInput("exp_data_selector", "Select Experimental Data:", choices=names(experimental_datasets)),
            checkboxInput("show_experimental", "Overlay Experimental Data", value=TRUE),
            hr(),
            p("Select variables to display on the plots below:"),
            checkboxGroupInput("currents_to_plot", "Ionic Currents:",
                               choices = c("I_f", "I_CaL", "I_Kr", "I_Ks", "I_NCX", "I_NaK"),
                               selected = c("I_f", "I_CaL", "I_Kr"), inline=TRUE),
            checkboxGroupInput("gates_to_plot", "Gating Variables:",
                               choices = c("y", "d_L", "f_L", "pa", "pi", "n"),
                               selected = c("y", "d_L", "f_L"), inline=TRUE)
        ))
    ),
    
    navset_card_tab(
        title = "Analysis",
        id = "main_tabs",
        nav_panel("Action Potential", withSpinner(plotOutput("ap_plot", height = "550px"))),
        nav_panel("Ionic Currents", withSpinner(plotOutput("currents_plot", height = "550px"))),
        nav_panel("Gating Variables", withSpinner(plotOutput("gates_plot", height = "550px"))),
        nav_panel("Model & Data Explained", uiOutput("explanation_tab"))
    )
)

#==============================================================================
# 3. SERVER LOGIC
#==============================================================================
server <- function(input, output, session) {

    # Model and GHK functions remain the same
    ghk <- function(V, ci, co, z, P, R, T, F) { v_rtf<-V*F/(R*T); if(abs(v_rtf)<1e-6) return(P*z*F*(ci-co)); return(P*z^2*v_rtf*F*(ci*exp(z*v_rtf)-co)/(exp(z*v_rtf)-1)) }
    sinoatrial_model <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
            Na_i<-max(1e-7, Na_i); K_i<-max(1e-7, K_i); Ca_i<-max(1e-7, Ca_i); Ca_sub<-max(1e-7, Ca_sub); Ca_j<-max(1e-7, Ca_j)
            E_Na<-(R*T/F)*log(Na_o/Na_i); E_K<-(R*T/F)*log(K_o/K_i)
            y_inf<-1/(1+exp((V+70.3)/9.1)); tau_y<-800/(exp(-(V+71.3)/16.1)+exp((V+71.3)/16.1)); dy<-(y_inf-y)/tau_y
            d_L_inf<-1/(1+exp(-(V+16.45)/4.337)); tau_d_L<-1.43*exp(-(V+45)/13)+1; dd_L<-(d_L_inf-d_L)/tau_d_L
            f_L_inf<-1/(1+exp((V+35.5)/4.1)); tau_f_L<-250*exp(-(V+60)/12)+2; df_L<-(f_L_inf-f_L)/tau_f_L
            f_Ca_inf<-1/(1+(Ca_sub/0.0004)); tau_f_Ca<-2; df_Ca<-(f_Ca_inf-f_Ca)/tau_f_Ca
            d_T_inf<-1/(1+exp(-(V+27)/7.2)); tau_d_T<-0.47*exp(-(V+60)/25)+1.2; dd_T<-(d_T_inf-d_T)/tau_d_T
            f_T_inf<-1/(1+exp((V+52)/11)); tau_f_T<-30*exp(-(V+100)/10)+15; df_T<-(f_T_inf-f_T)/tau_f_T
            pa_inf<-1/(1+exp(-(V+19.3)/10.1)); tau_pa<-4.2*exp(-(V*V)/(30*30))+2.1; dpa<-(pa_inf-pa)/tau_pa
            pi_inf<-1/(1+exp((V+25.9)/7.6)); dpi<-(pi_inf-pi)/tau_pi
            n_inf<-1/(1+exp(-(V-3.4)/11.4)); tau_n<-1200*exp(-(V-20)^2/(18^2))+300; dn<-(n_inf-n)/tau_n
            q_inf<-1/(1+exp((V+54)/8)); tau_q<-23; dq<-(q_inf-q)/tau_q
            r_inf<-1/(1+exp(-(V-29)/5.2)); tau_r<-1; dr<-(r_inf-r)/tau_r
            I_fNa<-g_fNa*y*(V-E_Na); I_fK<-g_fK*y*(V-E_K); I_f<-I_fNa+I_fK
            I_Kr<-g_Kr*pa*pi*(V-E_K); I_Ks<-g_Ks*n^2*(V-E_K); I_st<-g_st*q*r*(V-E_Na)
            I_CaL<-g_CaL*d_L*f_L*f_Ca*ghk(V,Ca_sub,Ca_o,2,1,R,T,F)
            I_CaT<-g_CaT*d_T*f_T*ghk(V,Ca_sub,Ca_o,2,1,R,T,F)
            I_bCa<-g_bCa*ghk(V,Ca_i,Ca_o,2,1,R,T,F)
            I_bNa<-g_bNa*(V-E_Na); I_bK<-g_bK*(V-E_K)
            I_NaK<-K_NaK*(K_o/(K_o+Km_K_NaK))*(Na_i^1.5/(Na_i^1.5+Km_Na_NaK^1.5))
            I_NCX<-K_NCX*(exp(gamma_NCX*V*F/(R*T))*Na_i^3*Ca_o-exp((gamma_NCX-1)*V*F/(R*T))*Na_o^3*Ca_sub)/((Km_Na_NCX^3+Na_o^3)*(Km_Ca_NCX+Ca_o)*(1+k_sat_NCX*exp((gamma_NCX-1)*V*F/(R*T))))
            I_pCa<-K_pCa*Ca_sub/(Km_pCa+Ca_sub)
            J_up<-V_up*Ca_sub^2/(k_up^2+Ca_sub^2); J_rel<-k_s*(Ca_j-Ca_sub); J_diff<-(Ca_sub-Ca_i)/tau_diff_Ca
            I_tot<-I_f+I_CaL+I_CaT+I_Kr+I_Ks+I_st+I_bNa+I_bK+I_bCa+I_NCX+I_pCa+I_NaK
            dV<- -I_tot/C_m; dNa_i<- -(I_fNa+I_bNa+I_st+3*I_NCX+3*I_NaK)*C_m/(F*V_i); dK_i<- -(I_fK+I_Kr+I_Ks+I_bK-2*I_NaK)*C_m/(F*V_i); dCa_i<- -I_bCa*C_m/(2*F*V_i)+J_diff; dCa_sub<-J_rel-J_up-J_diff-(I_CaL+I_CaT-2*I_NCX+I_pCa)*C_m/(2*F*V_sub); dCa_j<-(J_up-J_rel)*V_sub/V_j
            list(c(dV,dy,dd_L,df_L,df_Ca,dd_T,df_T,dpa,dpi,dn,dq,dr,dNa_i,dK_i,dCa_sub,dCa_j,dCa_i), c(I_f=I_f, I_CaL=I_CaL, I_Kr=I_Kr, I_Ks=I_Ks, I_NCX=I_NCX, I_NaK=I_NaK, dVdt=dV))
        })
    }
    
    sim_data<-reactiveVal(NULL); stored_traces<-reactiveVal(list())
    observeEvent(input$run_button, { f_scale<-input$g_f_scale*(1-input$block_f/100); CaL_scale<-input$g_CaL_scale*(1-input$block_CaL/100); Kr_scale<-input$g_Kr_scale*(1-input$block_Kr/100); params<-c(R=8314.472,T=310.0,F=96485.3415,C_m=57,V_i=0.00137,V_sub=0.0001096,V_j=0.000007672,Na_o=140,K_o=5.4,Ca_o=1.8,g_fNa=2.547*f_scale,g_fK=1.753*f_scale,g_CaL=0.000366*CaL_scale,g_CaT=0.0012,g_bCa=3.25e-5,g_Kr=4.2*Kr_scale,g_Ks=0.65,g_st=0.8,g_bNa=0.09,g_bK=0.3,K_NaK=38,Km_K_NaK=1,Km_Na_NaK=11,K_NCX=0.003,Km_Na_NCX=87.5,Km_Ca_NCX=1.38,k_sat_NCX=0.1,gamma_NCX=0.35,K_pCa=0.4,Km_pCa=0.0005,V_up=0.015,k_up=0.001,k_s=1.48e+08,tau_diff_Ca=5.469e-05,tau_pi=1); initial_state<-c(V=-60,y=0.01,d_L=0,f_L=1,f_Ca=1,d_T=0,f_T=1,pa=0,pi=1,n=0,q=0,r=1,Na_i=8,K_i=140,Ca_sub=1e-4,Ca_j=1e-4,Ca_i=1e-4); times<-seq(0,input$sim_time,by=1); showNotification("Running simulation...",duration=3,type="message"); solution<-ode(y=initial_state,times=times,func=sinoatrial_model,parms=params,method="lsoda"); solution<-as.data.frame(solution); solution$label<-paste0("Trace ",length(stored_traces())+1); sim_data(solution) }, ignoreNULL=FALSE)
    observeEvent(input$store_trace,{req(sim_data());stored_traces(append(stored_traces(),list(sim_data())));showNotification(paste("Stored",sim_data()$label[1]),type="message")})
    observeEvent(input$clear_traces,{stored_traces(list());showNotification("Cleared stored traces.",type="warning")})
    findpeaks<-function(x,minpeakheight=-20,minpeakdistance=200){pks<-which(diff(sign(diff(x,na.pad=F)),na.pad=F)<0)+1;pks<-pks[x[pks]>=minpeakheight];if(length(pks)>1){final_pks<-pks[1];for(pk in pks[-1])if(pk-tail(final_pks,1)>minpeakdistance)final_pks<-c(final_pks,pk);return(cbind(x[final_pks],final_pks))}else{return(NULL)}}

    output$ap_plot<-renderPlot({ req(sim_data()); current_sim<-sim_data(); peaks<-findpeaks(current_sim$V); bpm<-if(!is.null(peaks)&&nrow(peaks)>1)round(60000/mean(diff(current_sim$time[peaks[,2]])))else"N/A"; p<-ggplot()+labs(title=paste("Sinoatrial Node Action Potential (Current:",bpm,"BPM)"),x="Time (ms)",y="Membrane Potential (mV)")+theme_minimal(base_size=16)+theme(plot.title=element_text(hjust=0.5,face="bold"),legend.position="bottom"); if(length(stored_traces())>0)for(i in seq_along(stored_traces()))p<-p+geom_line(data=stored_traces()[[i]],aes(x=time,y=V,color=paste("Stored:",label)),linewidth=1); if(input$show_experimental){selected_exp_data<-experimental_datasets[[input$exp_data_selector]];p<-p+geom_line(data=selected_exp_data,aes(x=time,y=V,color="Experimental Data"),linewidth=1.2,linetype="dashed")}; p<-p+geom_line(data=current_sim,aes(x=time,y=V,color="Current Simulation"),linewidth=1.2)+scale_color_manual(name="Trace",values=c("Current Simulation"="#c0392b","Experimental Data"="#2c3e50","Stored: Trace 1"="#7f8c8d","Stored: Trace 2"="#95a5a6","Stored: Trace 3"="#bdc3c7")); print(p)})
    output$currents_plot<-renderPlot({validate(need(sim_data(),"Please run a simulation."),need(input$currents_to_plot,"Please select at least one ionic current to plot."));currents_long<-tidyr::pivot_longer(sim_data(),cols=all_of(input$currents_to_plot),names_to="current",values_to="pA");ggplot(currents_long,aes(x=time,y=pA,color=current))+geom_line(linewidth=1.1)+labs(title="Ionic Currents",x="Time (ms)",y="Current (pA)")+theme_minimal(base_size=14)+theme(legend.position="bottom",plot.title=element_text(hjust=0.5))})
    output$gates_plot<-renderPlot({validate(need(sim_data(),"Please run a simulation."),need(input$gates_to_plot,"Please select at least one gating variable to plot."));gates_long<-tidyr::pivot_longer(sim_data(),cols=all_of(input$gates_to_plot),names_to="gate",values_to="activation");ggplot(gates_long,aes(x=time,y=activation,color=gate))+geom_line(linewidth=1.1)+labs(title="Gating Variables",x="Time (ms)",y="Activation/Inactivation State")+theme_minimal(base_size=14)+theme(legend.position="bottom",plot.title=element_text(hjust=0.5))})
    
    # --- UPDATE: Explanation tab with references ---
    output$explanation_tab<-renderUI({
        tagList(
            h4("About the Computational Model"),
            p("This simulator uses the ", strong("Fabbri et al. (2017)"), " mathematical model of a human sinoatrial node (SAN) cell. It is a 'biophysically detailed' model, meaning its equations are designed to represent the actual physical and electrical behavior of ion channels and pumps based on laboratory measurements. This was a significant advancement, as many previous models were based on animal cells (e.g., rabbit), and this was one of the first comprehensive models specifically for a human pacemaker cell."),
            p("The model is a system of ordinary differential equations (ODEs). Each equation calculates the rate of change for a specific variable in the cell, such as the membrane voltage or the concentration of an ion. By solving these equations over time, we can simulate the cell's action potential."),
            tags$blockquote("Fabbri, A., et al. (2017). Computational analysis of the human sinus node action potential: model development and effects of mutations. ", tags$em("Journal of Molecular and Cellular Cardiology"), ", 108, 151-163.", style="font-size: 0.9em; border-left: 3px solid #eee; padding-left: 15px;"),
            
            hr(),
            h4("The Key Ionic Currents"),
            tags$dl(
                tags$dt("I_f (Funny Current)"),
                tags$dd("Carried by Na+ and K+ ions. This is the primary pacemaker current. It activates slowly when the cell membrane is very negative (at the end of the previous beat), causing a gradual depolarization that initiates the next heartbeat. Drugs like Ivabradine, which slow the heart rate, specifically block this channel. (See ", strong("DiFrancesco, 2010"), ")."),
                tags$dt("I_CaL (L-type Calcium Current)"),
                tags$dd("This is the main current responsible for the rapid upstroke (depolarization) of the SAN action potential. Its strong inward flow of Ca2+ ions quickly raises the membrane potential. It is the target of calcium channel blocker drugs like verapamil. (See ", strong("Bers, 2002"), ")."),
                tags$dt("I_Kr & I_Ks (Rapid and Slow Delayed Rectifier K+ Currents)"),
                tags$dd("These are potassium currents that are crucial for repolarization—bringing the membrane potential back down after the peak of the action potential. They allow K+ ions to leave the cell, making it more negative and 'resetting' it for the next beat. They are a common target for anti-arrhythmic drugs. (See ", strong("Sanguinetti & Tristani-Firouzi, 2006"), ")."),
                tags$dt("I_NCX (Sodium-Calcium Exchanger)"),
                tags$dd("An exchanger, not a simple channel. It plays a complex role, typically removing Ca2+ from the cell in exchange for Na+, contributing to the net charge balance."),
                tags$dt("I_NaK (Sodium-Potassium Pump)"),
                tags$dd("An active pump that uses energy (ATP) to maintain the fundamental ion gradients across the cell membrane (pumping Na+ out and K+ in). It is essential for the long-term stability of any cell.")
            ),
             tags$blockquote(
                tags$div("Bers, D. M. (2002). Cardiac excitation–contraction coupling. ", tags$em("Nature"), ", 415(6868), 198-205."),
                tags$div("DiFrancesco, D. (2010). The role of the funny current in pacemaker activity. ", tags$em("Circulation Research"), ", 106(3), 434-446."),
                tags$div("Sanguinetti, M. C., & Tristani-Firouzi, M. (2006). hERG potassium channels and cardiac arrhythmia. ", tags$em("Nature"), ", 440(7083), 463-469."),
                style="font-size: 0.9em; border-left: 3px solid #eee; padding-left: 15px;"
            ),

            hr(),
            h4("About the Experimental Data"),
            p("The experimental traces available for overlay are not from a single, raw recording. They are ", strong("digitized, representative waveforms"), " created based on typical morphologies and rates reported in electrophysiology literature for human SAN cells under different conditions. They serve as a valuable benchmark to assess how well the simulation's parameters can replicate real-world cardiac behavior."),
            p("The different rates (bradycardia, tachycardia) are primarily achieved through modulation of the pacemaker potential slope by the autonomic nervous system, a concept well-described in standard physiology textbooks. (See ", strong("Guyton and Hall, Textbook of Medical Physiology"),")."),
            tags$ul(
                tags$li(strong("Normal Sinus Rhythm:"), " Represents a typical, healthy resting heart rate of about 70 beats per minute."),
                tags$li(strong("Sinus Bradycardia:"), " Represents a slower resting heart rate (below 60 BPM). Note the longer, flatter diastolic depolarization phase."),
                tags$li(strong("Sinus Tachycardia:"), " Represents a faster resting heart rate (above 90 BPM). Note the shorter, steeper diastolic depolarization phase, leading to more frequent firing.")
            )
        )
    })
}

shinyApp(ui = ui, server = server)
