// This is a first attempt at writing an example of  a  mechanistic,
// compartmental, two-pool model in Rust language, utilizing the
// pharmsol crate to solve the ode's.  Created by JAM on Jan 11, 2026

// https://github.com/LAPKB/pharmsol

// Last updated on 13/1/26

// List necessary crates
use pharmsol::*;
use std::fs;


// MAIN function
fn main() {

    // First step is to build an experimental subject (model person)
    // An infusion bolus at time 0 sets the initial state:
    // For my two pool model, Pool  Sizes A: 6.0 B: 9.0, Thererfore Total: 15.0
    // An infusion is applied to pool A from time 0 to 10 with rate 7.0 (the FOA value).

    let subject = Subject::builder("1")
    // Apply infusion from time 0.0 to 10.0 [0,10]  into Pool A (compartment 0) at a constant  rate of 7.0 over [0, 10]
	.infusion(0.0, 10.0, 0, 7.0)
    //Observations are scheduled at selected timepoints.
        .observation(0.0, 0.0, 0)
        .repeat(19, 1.0)
        .observation(0.0, 0.0, 1)
        .repeat(19, 1.0)
        .observation(0.0, 0.0, 2)
        .repeat(19, 1.0)
        .build();

    // Parameter values corresponding to:
    // SA = 20.0, SB = 25.0, VAB = 18.0, VBA = 13.0, VBO = 8.0, KAB = 0.32, KBA = 0.36, KBO = 0.31
    let params = vec![20.0, 25.0, 18.0, 13.0, 8.0, 0.32, 0.36, 0.31];

    // Create the ODE model for the two-pool system.
    // The  state vector x has three components: x[0]=QA, x[1]=QB, x[2]=QT (Total).
    // The infusion rate is passed in via rateiv[0] (applied to pool A).

    let ode = equation::ODE::new(
        // ODE function:
        |x, p, _t, dx, rateiv, _cov| {

            // Extract parameters (order: SA, SB, VAB, VBA, VBO, KAB, KBA, KBO)
            fetch_params!(p, sa, sb, vab, vba, vbo, kab, kba, kbo);

            // Retrieve state values
            let qa = x[0];
            let qb = x[1];

            // Compute the concentrations in each pool
            let con_a = qa / sa;
            let con_b = qb / sb;

            // Calculate the HMM fluxes
            let fab = vab / (1.0 + (kab / con_a));
            let fba = vba / (1.0 + (kba / con_b));
            let fbo = vbo / (1.0 + (kbo / con_b));

            // Differential equations, using the infusion rate provided as rateiv[0]
            dx[0] = rateiv[0] + fba - fab; // pool A
            dx[1] = fab - fba - fbo; // pool B
            dx[2] = rateiv[0] - fbo; // total pool
        },

	// The remaining closures are for lag, analytic approximation and covariance handling.
        |_p| lag! {},
        |_p| fa! {},
        |_p, _t, _cov, x| {
            x[0] = 6.0; // Initial QA
            x[1] = 9.0; // Initial QB	    
            x[2] = 15.0; // Initial QT
        },

        // Observation function: here we output the state vector directly.
        |x, _p, _t, _cov, y| {
            y[0] = x[0];
            y[1] = x[1];
            y[2] = x[2];
        },


	// create function to output  data to file


	
	// specify dimensions: three states and three observations
        (3, 3), 
    );

      // Estimate predictions from the ODE model using the subject data and parameters.
    let ode_predictions = ode.estimate_predictions(&subject, &params);

    // Assuming that `ode_predictions.get_predictions()` returns a vector of Prediction structs.
    let predictions = ode_predictions.get_predictions();

    // Create a HashMap to group predictions by the outeq value.
    let mut groups: HashMap<usize, Vec<f64>> = HashMap::new();

    for pred in predictions.iter() {
        // Insert the prediction into the group corresponding to its outeq.
        groups
            .entry(pred.outeq())
            .or_insert_with(Vec::new)
            .push(pred.prediction());
    }

    // Print the grouped predictions to screen.
    for (outeq, preds) in &groups {
        println!("Output Equation {} has {} predictions:", outeq, preds.len());
        for pred in preds {
            println!("{:#?}", pred);
        }
	
    }
    //output predictions to a file
    // for (outeq, preds) in &groups {
    // 	//	fs::write(outeq, "./output.dat");
    // 	fs::write("./output.dat", outeq); 
    // 	}
    
    }
