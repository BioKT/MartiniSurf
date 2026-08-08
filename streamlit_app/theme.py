from __future__ import annotations

import streamlit as st


def apply_theme() -> None:
    st.markdown(
        """
        <style>
        :root {
          --ms-ink: #17211c;
          --ms-muted: #63756d;
          --ms-line: #ecd4df;
          --ms-panel: #fff8fb;
          --ms-accent: #ff4fa7;
          --ms-teal: #2d9bb3;
          --ms-deep: #16312f;
        }
        .stApp {
          background:
            radial-gradient(circle at top left, rgba(255, 79, 167, 0.12), transparent 30rem),
            linear-gradient(180deg, #fffafd 0%, #eef8f8 100%);
          color: var(--ms-ink);
        }
        [data-testid="stSidebar"] {
          background: #16312f;
        }
        [data-testid="stSidebar"] * {
          color: #f5fbf7;
        }
        h1, h2, h3 {
          letter-spacing: 0;
        }
        h1 {
          font-size: clamp(3rem, 6vw, 5.8rem) !important;
          line-height: 0.92 !important;
        }
        h2 {
          color: var(--ms-deep);
        }
        div[data-testid="stMetric"] {
          background: rgba(255,255,255,0.72);
          border: 1px solid var(--ms-line);
          border-radius: 8px;
          padding: 0.75rem 0.9rem;
        }
        .ms-header {
          border-bottom: 1px solid var(--ms-line);
          padding: 0.4rem 0 1.35rem 0;
          margin-bottom: 1.1rem;
          display: grid;
          grid-template-columns: minmax(0, 1fr) 7.5rem;
          gap: 1.25rem;
          align-items: center;
        }
        .ms-kicker {
          color: var(--ms-teal);
          font-size: 0.82rem;
          font-weight: 700;
          text-transform: uppercase;
        }
        .ms-brand {
          color: var(--ms-accent);
          font-size: clamp(3.6rem, 8vw, 7.4rem);
          font-weight: 900;
          line-height: 0.88;
          margin: 0.18rem 0 0.45rem;
        }
        .ms-soft {
          color: var(--ms-muted);
        }
        .ms-hero-logo img {
          max-height: 12rem;
          object-fit: contain;
        }
        .ms-chip-row {
          display: flex;
          gap: 0.45rem;
          flex-wrap: wrap;
          margin-top: 0.8rem;
        }
        .ms-chip {
          border: 1px solid var(--ms-line);
          border-radius: 999px;
          color: var(--ms-deep);
          background: rgba(255,255,255,0.72);
          padding: 0.28rem 0.68rem;
          font-size: 0.82rem;
          font-weight: 650;
        }
        .ms-section {
          border: 1px solid var(--ms-line);
          border-radius: 8px;
          background: rgba(255,255,255,0.62);
          padding: 1rem;
          margin-bottom: 1rem;
        }
        .ms-log-row {
          border-left: 3px solid var(--ms-accent);
          padding: 0.35rem 0 0.35rem 0.75rem;
          margin: 0.35rem 0;
          background: rgba(255,255,255,0.45);
        }
        .ms-log-label {
          color: var(--ms-teal);
          font-size: 0.78rem;
          font-weight: 750;
          text-transform: uppercase;
        }
        .ms-log-value {
          color: var(--ms-ink);
          font-size: 0.92rem;
          margin-top: 0.08rem;
        }
        .stButton > button, .stDownloadButton > button {
          border-radius: 8px;
          border: 1px solid var(--ms-accent);
        }
        .stButton > button[kind="primary"] {
          background: var(--ms-accent);
        }
        div[data-testid="stTabs"] button[aria-selected="true"] {
          color: var(--ms-accent);
        }
        @media (max-width: 720px) {
          .ms-header {
            grid-template-columns: 1fr;
          }
          .ms-hero-logo {
            display: none;
          }
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
