use eframe::egui;
use nalgebra::Vector2;
use rand::Rng;
use rayon::prelude::*;
use std::time::Instant;

const G: f64 = 6.67430e-3;

#[derive(Clone)]
struct Body {
    pos: Vector2<f64>,
    vel: Vector2<f64>,
    mass: f64,
}

impl Body {
    fn update(&mut self, f: Vector2<f64>, dt: f64) {
        self.vel += f / self.mass * dt;
        self.pos += self.vel * dt;
    }
}

fn compute_force(a: &Body, bodies: &[Body]) -> Vector2<f64> {
    let mut f = Vector2::zeros();
    for b in bodies {
        if (a.pos - b.pos).magnitude() > 1e-5 {
            let dir = b.pos - a.pos;
            let d = dir.magnitude() + 1e-3;
            f += dir.normalize() * (G * a.mass * b.mass / (d * d));
        }
    }
    f
}

struct Sim {
    bodies: Vec<Body>,
    last_update: Instant,
    paused: bool,
    dt: f64,
    scale: f64,
}

impl Default for Sim {
    fn default() -> Self {
        let mut rng = rand::thread_rng();
        let bodies = (0..300)
            .map(|_| Body {
                pos: Vector2::new(
                    rng.gen_range(-400.0..400.0),
                    rng.gen_range(-400.0..400.0),
                ),
                vel: Vector2::new(
                    rng.gen_range(-2.0..2.0),
                    rng.gen_range(-2.0..2.0),
                ),
                mass: rng.gen_range(1.0..10.0),
            })
            .collect();

        Self {
            bodies,
            last_update: Instant::now(),
            paused: false,
            dt: 0.1,
            scale: 0.5,
        }
    }
}

impl eframe::App for Sim {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::TopBottomPanel::top("controls").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if ui.button(if self.paused { "▶️ Resume" } else { "⏸ Pause" }).clicked() {
                    self.paused = !self.paused;
                }
                ui.label("Δt:");
                ui.add(egui::Slider::new(&mut self.dt, 0.001..=1.0).logarithmic(true));
                ui.label("Zoom:");
                ui.add(egui::Slider::new(&mut self.scale, 0.1..=1.0));
            });
        });

        let now = Instant::now();
        let elapsed = now.duration_since(self.last_update).as_secs_f64();
        self.last_update = now;

        if !self.paused {
            // Compute forces in parallel
            let forces: Vec<_> = self
                .bodies
                .par_iter()
                .map(|b| compute_force(b, &self.bodies))
                .collect();

            // Update bodies
            self.bodies
                .par_iter_mut()
                .zip(forces)
                .for_each(|(b, f)| b.update(f, self.dt * elapsed * 60.0)); // scale for frame rate
        }

        // Draw all particles
        egui::CentralPanel::default().show(ctx, |ui| {
            let (response, painter) =
                ui.allocate_painter(ui.available_size(), egui::Sense::hover());

            let center = response.rect.center();
            let scale = self.scale as f32;

            for b in &self.bodies {
                let pos = egui::pos2(
                    center.x + (b.pos.x as f32 * scale),
                    center.y + (b.pos.y as f32 * scale),
                );
                painter.circle_filled(pos, 1.5, egui::Color32::WHITE);
            }
        });

        ctx.request_repaint(); // ensure continuous frames
    }
}

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([900.0, 900.0])
            .with_title("N-Body Simulation (Egui)"),
        ..Default::default()
    };

    eframe::run_native(
        "N-Body Simulation",
        options,
        Box::new(|_cc| Box::<Sim>::default()),
    )
}
