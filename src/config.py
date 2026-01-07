from dataclasses import dataclass

@dataclass(frozen=True)
class Config:
    random_state: int = 42
    test_size: float = 0.2
