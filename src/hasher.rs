use std::hash::{BuildHasher, Hasher};
use xxhash_rust::xxh3::Xxh3;

#[derive(Clone)]
pub struct Xxh3Builder {
    pub seed: u64,
}

impl BuildHasher for Xxh3Builder {
    type Hasher = Xxh3Hasher;

    fn build_hasher(&self) -> Self::Hasher {
        Xxh3Hasher {
            inner: Xxh3::with_seed(self.seed),
        }
    }
}

pub struct Xxh3Hasher {
    inner: Xxh3,
}

impl Hasher for Xxh3Hasher {
    fn write(&mut self, bytes: &[u8]) {
        self.inner.update(bytes);
    }

    fn finish(&self) -> u64 {
        self.inner.clone().digest()
    }
}
