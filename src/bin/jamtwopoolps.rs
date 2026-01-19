// This is a first attempt at writing an example of  a  mechanistic,
// compartmental, two-pool model in Rust language, utilizing the
// pharmsol crate to solve the ode's.

//Created by JAM on Jan 11, 2026 at Norwich UK

// https://github.com/LAPKB/pharmsol

// Last updated on 16/1/26

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

    // Create CSV output: each row is a time point, columns are time and each output equation
    // Collect all predictions with their time and outeq
    let mut all_preds: Vec<(f64, usize, f64)> = predictions
        .iter()
        .map(|p| (p.time(), p.outeq(), p.prediction()))
        .collect();

    // Sort by time, then by outeq
    all_preds.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.1.cmp(&b.1))
    });

    // Get unique sorted time points
    let mut times: Vec<f64> = all_preds.iter().map(|(t, _, _)| *t).collect();
    times.dedup();

    // Determine the number of output equations
    let num_outeqs = groups.keys().max().map_or(0, |&m| m + 1);

    // Build CSV content
    let mut csv_content = String::new();

    // Header row
    csv_content.push_str("Time,QA,QB,QT\n");

    // Data rows - group predictions by time
    for time in &times {
        csv_content.push_str(&format!("{:.4}", time));

        // Find all predictions at this time point
        let preds_at_time: Vec<_> = all_preds
            .iter()
            .filter(|(t, _, _)| (*t - time).abs() < 1e-10)
            .collect();

        // Create row with predictions for each outeq
        for outeq in 0..num_outeqs {
            if let Some((_, _, pred_val)) = preds_at_time.iter().find(|(_, o, _)| *o == outeq) {
                csv_content.push_str(&format!(",{:.6}", pred_val));
            } else {
                csv_content.push_str(",");
            }
        }
        csv_content.push_str("\n");
    }

    // Write to file
    fs::write("./predictions.csv", &csv_content).expect("Failed to write CSV file");
    println!("\nPredictions saved to predictions.csv");

    // Also print to console
    println!("\nCSV Output:\n{}", csv_content);
}
