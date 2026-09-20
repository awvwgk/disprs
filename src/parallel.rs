use std::panic::resume_unwind;
use std::sync::OnceLock;

/// Resolved once per process; `DISPRS_NUM_THREADS` is not re-read after the first call.
fn default_threads() -> usize {
    static THREADS: OnceLock<usize> = OnceLock::new();
    *THREADS.get_or_init(|| {
        std::env::var("DISPRS_NUM_THREADS")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or_else(|| {
                std::thread::available_parallelism()
                    .map_or(1, usize::from)
                    .min(4)
            })
    })
}

pub(crate) fn thread_count(work: usize, minimum: usize) -> usize {
    if work < minimum {
        return 1;
    }
    default_threads().clamp(1, work.max(1))
}

pub(crate) fn map<T: Send>(
    work: usize,
    minimum: usize,
    operation: impl Fn(usize, usize) -> T + Sync,
) -> Vec<T> {
    map_with_threads(thread_count(work, minimum), operation)
}

fn map_with_threads<T: Send>(
    threads: usize,
    operation: impl Fn(usize, usize) -> T + Sync,
) -> Vec<T> {
    if threads == 1 {
        return vec![operation(0, 1)];
    }
    std::thread::scope(|scope| {
        let operation = &operation;
        let handles: Vec<_> = (1..threads)
            .map(|worker| scope.spawn(move || operation(worker, threads)))
            .collect();
        let mut results = vec![operation(0, threads)];
        results.extend(handles.into_iter().map(|handle| {
            handle
                .join()
                .unwrap_or_else(|payload| resume_unwind(payload))
        }));
        results
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cyclic_rows_cover_each_pair_and_triplet_once() {
        for atoms in [0, 1, 2, 7, 64] {
            for threads in [1, 2, 4] {
                let rows = map_with_threads(threads, |start, stride| {
                    (start..atoms).step_by(stride).collect::<Vec<_>>()
                });
                let mut rows: Vec<_> = rows.into_iter().flatten().collect();
                rows.sort_unstable();
                assert_eq!(rows, (0..atoms).collect::<Vec<_>>());
            }
        }
    }
}
