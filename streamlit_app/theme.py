from __future__ import annotations

import streamlit as st


def apply_theme() -> None:
    st.markdown(
        """
        <style>
        :root {
          --ms-bg: #0E0D11;
          --ms-sidebar-bg: #121016;
          --ms-bg-elevated: #19161E;
          --ms-bg-panel: #211C27;
          --ms-bg-input: #17141B;
          --ms-bg-hover: #2B2432;
          --ms-border: rgba(229, 214, 226, 0.16);
          --ms-border-strong: rgba(255, 108, 183, 0.50);
          --ms-text: #FAF5F8;
          --ms-text-secondary: #D6CAD3;
          --ms-text-muted: #A99DA7;
          --ms-text-on-accent: #180710;
          --ms-primary: #FF62B4;
          --ms-primary-hover: #FF9ED2;
          --ms-success: #68D391;
          --ms-warning: #F2B86B;
          --ms-error: #F07883;
          --ms-focus: rgba(255, 98, 180, 0.34);
          --ms-accent: #FF62B4;
          --ms-radius-sm: 8px;
          --ms-radius-md: 12px;
          --ms-radius-lg: 16px;
          --ms-shadow-soft: 0 14px 36px rgba(0, 0, 0, 0.24);
          --ms-shadow-panel: 0 10px 28px rgba(0, 0, 0, 0.18);

          --ms-ink: var(--ms-text);
          --ms-muted: var(--ms-text-secondary);
          --ms-line: var(--ms-border);
          --ms-panel: rgba(33, 28, 39, 0.92);
          --ms-teal: #8CDADF;
          --ms-mint: #A5E8DC;
          --ms-deep: #0E0D11;
          --ms-green: var(--ms-success);
          --ms-amber: var(--ms-warning);
        }
        *,
        *::before,
        *::after {
          box-sizing: border-box;
        }
        .stApp {
          background:
            linear-gradient(180deg, rgba(255, 98, 180, 0.045), transparent 18rem),
            linear-gradient(145deg, #0E0D11 0%, #141118 52%, #0A090D 100%);
          color: var(--ms-text);
          font-family: Inter, "Source Sans Pro", system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
        }
        [data-testid="stAppViewContainer"],
        [data-testid="stMain"] {
          background: transparent !important;
        }
        html,
        body,
        #root,
        .stApp,
        [data-testid="stAppViewContainer"] {
          min-height: 100vh !important;
          overflow-y: auto !important;
        }
        .stApp > header,
        [data-testid="stHeader"],
        header[data-testid="stHeader"] {
          display: block !important;
          height: 3rem !important;
          min-height: 3rem !important;
          opacity: 1 !important;
          pointer-events: auto !important;
          background: transparent !important;
          background-color: transparent !important;
          box-shadow: none !important;
          overflow: visible !important;
          z-index: 999998 !important;
        }
        [data-testid="stHeader"]::before,
        [data-testid="stHeader"]::after {
          display: none !important;
          background: transparent !important;
        }
        [data-testid="stDecoration"],
        [data-testid="stStatusWidget"] {
          display: none !important;
        }
        [data-testid="stToolbar"] {
          background: transparent !important;
          pointer-events: auto !important;
        }
        [data-testid="collapsedControl"],
        [data-testid="stSidebarCollapsedControl"],
        button[data-testid="stBaseButton-headerNoPadding"],
        button[title="Open sidebar"],
        button[aria-label="Open sidebar"],
        button[title="Close sidebar"],
        button[aria-label="Close sidebar"] {
          display: flex !important;
          visibility: visible !important;
          opacity: 1 !important;
          pointer-events: auto !important;
          background: rgba(25, 22, 30, 0.96) !important;
          border: 1px solid rgba(255, 98, 180, 0.48) !important;
          border-radius: 8px !important;
          color: var(--ms-ink) !important;
          z-index: 999999 !important;
        }
        [data-testid="collapsedControl"] *,
        [data-testid="stSidebarCollapsedControl"] *,
        button[data-testid="stBaseButton-headerNoPadding"] *,
        button[title="Open sidebar"] *,
        button[aria-label="Open sidebar"] * {
          color: #ffffff !important;
          fill: #ffffff !important;
          stroke: #ffffff !important;
        }
        [data-testid="collapsedControl"] svg,
        [data-testid="stSidebarCollapsedControl"] svg,
        button[data-testid="stBaseButton-headerNoPadding"] svg,
        button[title="Open sidebar"] svg,
        button[aria-label="Open sidebar"] svg,
        button[title="Close sidebar"] svg,
        button[aria-label="Close sidebar"] svg {
          color: #ffffff !important;
          fill: #ffffff !important;
          stroke: #ffffff !important;
        }
        [data-testid="stSidebar"] {
          background: linear-gradient(180deg, #121016 0%, #0C0B10 100%);
          border-right: 1px solid var(--ms-line);
          scrollbar-color: rgba(255, 98, 180, 0.78) rgba(20, 17, 24, 0.75);
          scrollbar-width: thin;
        }
        [data-testid="stSidebarContent"] {
          overflow-y: auto !important;
          max-height: 100vh !important;
          padding-top: 1rem;
          padding-left: 1.05rem;
          padding-right: 1.05rem;
          scrollbar-color: rgba(255, 98, 180, 0.78) rgba(20, 17, 24, 0.75);
          scrollbar-width: thin;
        }
        [data-testid="stSidebar"]::-webkit-scrollbar,
        [data-testid="stSidebarContent"]::-webkit-scrollbar,
        [data-testid="stAppViewContainer"]::-webkit-scrollbar,
        body::-webkit-scrollbar {
          width: 10px;
          height: 10px;
        }
        [data-testid="stSidebar"]::-webkit-scrollbar-track,
        [data-testid="stSidebarContent"]::-webkit-scrollbar-track,
        [data-testid="stAppViewContainer"]::-webkit-scrollbar-track,
        body::-webkit-scrollbar-track {
          background: rgba(20, 17, 24, 0.75);
        }
        [data-testid="stSidebar"]::-webkit-scrollbar-thumb,
        [data-testid="stSidebarContent"]::-webkit-scrollbar-thumb,
        [data-testid="stAppViewContainer"]::-webkit-scrollbar-thumb,
        body::-webkit-scrollbar-thumb {
          background: linear-gradient(180deg, var(--ms-teal), var(--ms-accent));
          border: 2px solid rgba(20, 17, 24, 0.75);
          border-radius: 999px;
        }
        [data-testid="stSidebar"]::-webkit-scrollbar-thumb:hover,
        [data-testid="stSidebarContent"]::-webkit-scrollbar-thumb:hover,
        [data-testid="stAppViewContainer"]::-webkit-scrollbar-thumb:hover,
        body::-webkit-scrollbar-thumb:hover {
          background: var(--ms-mint, #64E3C4);
        }
        [data-testid="stSidebar"] {
          color: var(--ms-text);
        }
        [data-testid="stSidebar"] p,
        [data-testid="stSidebar"] label,
        [data-testid="stSidebar"] span {
          color: inherit;
        }
        [data-testid="stSidebar"] img {
          background: linear-gradient(180deg, rgba(33, 28, 39, 0.96), rgba(18, 16, 22, 0.90));
          border: 1px solid rgba(255, 98, 180, 0.22);
          border-radius: var(--ms-radius-lg);
          padding: 0.55rem;
          height: 13.5rem;
          object-fit: contain;
          box-shadow: var(--ms-shadow-panel);
        }
        .block-container {
          padding-top: 3rem;
          padding-bottom: 2rem;
          max-width: 1500px;
        }
        h1, h2, h3 {
          letter-spacing: 0;
        }
        h2 {
          color: var(--ms-ink);
        }
        div[data-testid="stMetric"] {
          background: linear-gradient(180deg, rgba(35, 30, 41, 0.94), rgba(24, 21, 29, 0.92));
          border: 1px solid var(--ms-line);
          border-radius: var(--ms-radius-md);
          padding: 0.58rem 0.78rem;
          min-height: 4.4rem;
          box-shadow: 0 6px 18px rgba(0,0,0,0.13);
        }
        div[data-testid="stMetric"] label,
        div[data-testid="stMetric"] [data-testid="stMetricValue"] {
          color: var(--ms-text);
        }
        div[data-testid="stMetric"] label {
          color: var(--ms-text-secondary) !important;
          font-size: 0.76rem !important;
          letter-spacing: 0;
        }
        div[data-testid="stMetric"] [data-testid="stMetricValue"] {
          font-family: "Roboto Mono", "SFMono-Regular", Consolas, monospace;
          font-size: 1.1rem !important;
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
        .ms-brand-small {
          color: var(--ms-primary);
          font-size: clamp(2rem, 4vw, 3.8rem);
          font-weight: 950;
          line-height: 1;
          margin: 0.25rem 0 0.35rem;
        }
        .ms-soft {
          color: var(--ms-muted);
        }
        .ms-sidebar-title {
          color: var(--ms-text);
          font-size: 1.22rem;
          font-weight: 900;
          line-height: 1;
          margin: 0;
          letter-spacing: 0;
        }
        .ms-sidebar-title span {
          color: var(--ms-accent);
        }
        .ms-sidebar-title strong {
          color: var(--ms-primary);
          font-weight: 900;
        }
        .ms-topbar {
          display: grid;
          grid-template-columns: minmax(18rem, 1fr) auto;
          gap: 1rem;
          align-items: center;
          border: 1px solid var(--ms-line);
          background: linear-gradient(180deg, rgba(33, 28, 39, 0.93), rgba(17, 15, 21, 0.90));
          border-radius: var(--ms-radius-lg);
          padding: 0.82rem 1rem;
          margin-bottom: 1.05rem;
          box-shadow: var(--ms-shadow-soft);
        }
        .ms-top-brand {
          display: flex;
          gap: 0;
          align-items: center;
          min-width: 0;
        }
        .ms-top-brand img {
          width: 3.8rem;
          height: 3.8rem;
          object-fit: contain;
          background: rgba(24, 21, 29, 0.86);
          border-radius: var(--ms-radius-md);
          padding: 0.28rem;
        }
        .ms-top-title {
          color: var(--ms-primary);
          font-size: clamp(1.65rem, 2.6vw, 2.65rem);
          font-weight: 950;
          line-height: 0.98;
          letter-spacing: 0;
        }
        .ms-top-subtitle {
          color: var(--ms-muted);
          font-size: 0.88rem;
          margin-top: 0.2rem;
          white-space: nowrap;
          overflow: hidden;
          text-overflow: ellipsis;
          max-width: 34rem;
        }
        .ms-top-pills {
          display: flex;
          gap: 0.5rem;
          flex-wrap: wrap;
          justify-content: flex-end;
        }
        .ms-top-pills span,
        .ms-ready {
          border: 1px solid var(--ms-line);
          border-radius: 999px;
          background: rgba(255,255,255,0.045);
          padding: 0.42rem 0.68rem;
          color: var(--ms-ink);
          font-size: 0.82rem;
          font-weight: 760;
        }
        .ms-top-pills span:nth-child(2) {
          color: var(--ms-teal);
          border-color: rgba(140, 218, 223, 0.40);
        }
        .ms-top-pills span:nth-child(3) {
          color: var(--ms-green);
          border-color: rgba(91, 234, 137, 0.35);
        }
        .ms-ready {
          display: flex;
          gap: 0.45rem;
          align-items: center;
          justify-content: center;
        }
        .ms-model-toggles + div[data-testid="stHorizontalBlock"] div[data-testid="stToggle"] {
          border: 1px solid var(--ms-line);
          background: rgba(33, 28, 39, 0.82);
          border-radius: var(--ms-radius-lg);
          min-height: 4.25rem;
          padding: 0.9rem 1rem;
        }
        .ms-model-toggles + div[data-testid="stHorizontalBlock"] div[data-testid="stToggle"] label p {
          font-size: 1.02rem;
          font-weight: 850;
        }
        .ms-ready span {
          width: 0.62rem;
          height: 0.62rem;
          border-radius: 50%;
          background: var(--ms-green);
          box-shadow: 0 0 12px rgba(85, 217, 138, 0.58);
        }
        .ms-steps {
          margin: 1.4rem 0 1rem;
        }
        .ms-step {
          display: flex;
          align-items: center;
          gap: 0.7rem;
          padding: 0.58rem 0.2rem;
          color: var(--ms-muted);
        }
        .ms-step span {
          width: 1.82rem;
          height: 1.82rem;
          border-radius: 50%;
          border: 1px solid rgba(255, 98, 180, 0.42);
          background: rgba(255,255,255,0.04);
        }
        .ms-step p {
          margin: 0;
          font-weight: 700;
        }
        .ms-step.active {
          color: var(--ms-deep);
          background: linear-gradient(135deg, var(--ms-primary), var(--ms-primary-hover));
          border-radius: var(--ms-radius-md);
          padding-left: 0.7rem;
          margin: 0.25rem 0;
        }
        .ms-step.active span {
          border-color: var(--ms-deep);
        }
        .ms-step.done span {
          background: var(--ms-green);
          border-color: var(--ms-green);
        }
        [data-testid="stSidebar"] div[role="radiogroup"] {
          display: grid;
          gap: 0.35rem;
        }
        [data-testid="stSidebar"] div[role="radiogroup"] label {
          border: 1px solid transparent;
          border-radius: var(--ms-radius-md);
          padding: 0.52rem 0.58rem;
          margin: 0;
          min-height: 2.65rem;
          transition: background-color 160ms ease, border-color 160ms ease, color 160ms ease;
        }
        [data-testid="stSidebar"] div[role="radiogroup"] label:hover {
          background: rgba(255, 98, 180, 0.09);
          border-color: rgba(255, 98, 180, 0.20);
        }
        [data-testid="stSidebar"] div[role="radiogroup"] label:has(input:checked) {
          background: linear-gradient(90deg, rgba(255, 92, 173, 0.18), rgba(126, 221, 226, 0.07));
          border-color: rgba(255, 98, 180, 0.46);
          color: var(--ms-text) !important;
          font-weight: 820;
        }
        .ms-control-title,
        .ms-panel-title {
          color: var(--ms-ink);
          font-size: 1.12rem;
          font-weight: 850;
          margin: 0.25rem 0 0.8rem;
          letter-spacing: 0;
        }
        div[data-testid="stTabs"] {
          border: 1px solid var(--ms-line);
          border-radius: var(--ms-radius-lg);
          background: rgba(27, 23, 32, 0.78);
          padding: 0.75rem;
        }
        div[data-testid="stTabs"] button {
          color: var(--ms-muted) !important;
        }
        div[data-testid="stTabs"] button p {
          color: inherit !important;
          font-weight: 750;
        }
        div[data-testid="stTabs"] button[aria-selected="true"] {
          color: var(--ms-deep) !important;
          background: linear-gradient(135deg, var(--ms-primary), var(--ms-primary-hover));
          border-radius: var(--ms-radius-md);
        }
        div[data-testid="stRadio"] label,
        div[data-testid="stCheckbox"] label,
        div[data-testid="stToggle"] label {
          color: var(--ms-text) !important;
        }
        div[data-testid="stRadio"] label:hover,
        div[data-testid="stCheckbox"] label:hover,
        div[data-testid="stToggle"] label:hover {
          color: var(--ms-primary-hover) !important;
        }
        div[data-testid="stProgress"] > div > div {
          background-color: var(--ms-primary) !important;
        }
        div[data-testid="stProgress"] p {
          color: #06111A !important;
          font-weight: 850 !important;
        }
        [data-testid="stWidgetLabel"] p,
        [data-testid="stWidgetLabel"] label,
        [data-testid="stWidgetLabel"] span {
          color: #D8E6EC !important;
          font-size: 0.83rem;
          font-weight: 760;
          letter-spacing: 0;
        }
        div[data-testid="stTextInput"] input,
        div[data-testid="stTextArea"] textarea,
        div[data-testid="stNumberInput"] input {
          background-color: var(--ms-bg-input) !important;
          color: var(--ms-text) !important;
          border: 1px solid var(--ms-border) !important;
          border-radius: var(--ms-radius-md) !important;
          box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.03) !important;
          min-height: 2.55rem;
        }
        div[data-testid="stTextInput"] input:hover,
        div[data-testid="stTextArea"] textarea:hover,
        div[data-testid="stNumberInput"] input:hover {
          border-color: var(--ms-border-strong) !important;
          background-color: #211B27 !important;
        }
        div[data-testid="stTextInput"] input:focus,
        div[data-testid="stTextArea"] textarea:focus,
        div[data-testid="stNumberInput"] input:focus {
          border-color: var(--ms-primary) !important;
          box-shadow: 0 0 0 0.18rem var(--ms-focus) !important;
        }
        div[data-baseweb="select"] > div,
        div[data-testid="stMultiSelect"] div[data-baseweb="select"] > div {
          background-color: var(--ms-bg-input) !important;
          color: var(--ms-text) !important;
          border-color: var(--ms-border) !important;
          border-radius: var(--ms-radius-md) !important;
          min-height: 2.55rem;
        }
        div[data-baseweb="select"] span,
        div[data-baseweb="select"] svg {
          color: var(--ms-text-secondary) !important;
          fill: var(--ms-text-secondary) !important;
        }
        div[data-baseweb="popover"],
        div[data-baseweb="menu"],
        ul[role="listbox"] {
          background: var(--ms-bg-elevated) !important;
          color: var(--ms-text) !important;
          border: 1px solid var(--ms-border-strong) !important;
          border-radius: var(--ms-radius-md) !important;
          box-shadow: 0 18px 50px rgba(0,0,0,0.42) !important;
        }
        li[role="option"],
        div[role="option"] {
          background: transparent !important;
          color: var(--ms-text) !important;
        }
        li[role="option"]:hover,
        div[role="option"]:hover,
        li[role="option"][aria-selected="true"],
        div[role="option"][aria-selected="true"] {
          background: #2B2432 !important;
          color: var(--ms-text) !important;
        }
        div[data-baseweb="tag"] {
          background: rgba(255, 98, 180, 0.16) !important;
          border: 1px solid rgba(255, 98, 180, 0.36) !important;
          color: var(--ms-text) !important;
        }
        div[data-baseweb="tag"] span,
        div[data-baseweb="tag"] svg {
          color: var(--ms-text) !important;
          fill: var(--ms-text) !important;
        }
        input::placeholder,
        textarea::placeholder {
          color: var(--ms-text-muted) !important;
        }
        [data-testid="stFileUploader"] {
          margin-bottom: 0.8rem;
        }
        [data-testid="stFileUploader"] section,
        [data-testid="stFileUploaderDropzone"] {
          background: linear-gradient(180deg, rgba(27, 23, 32, 0.96), rgba(16, 14, 20, 0.96)) !important;
          border: 1px dashed rgba(255, 98, 180, 0.34) !important;
          border-radius: var(--ms-radius-lg) !important;
          min-height: 4.55rem;
          box-shadow: inset 0 0 0 1px rgba(255,255,255,0.015);
        }
        [data-testid="stFileUploader"] section p,
        [data-testid="stFileUploader"] section span,
        [data-testid="stFileUploader"] section small {
          color: var(--ms-text-secondary) !important;
        }
        [data-testid="stFileUploader"] button {
          background: rgba(255, 98, 180, 0.12) !important;
          border: 1px solid rgba(255, 98, 180, 0.38) !important;
          color: var(--ms-ink) !important;
          border-radius: var(--ms-radius-md) !important;
        }
        [data-testid="stExpander"] {
          border: 1px solid var(--ms-line);
          background: rgba(33, 28, 39, 0.82);
          border-radius: var(--ms-radius-lg);
          box-shadow: 0 8px 22px rgba(0,0,0,0.12);
        }
        [data-testid="stExpander"] summary,
        [data-testid="stExpander"] summary p {
          color: var(--ms-text) !important;
          font-weight: 750;
        }
        [data-testid="stMarkdownContainer"] code,
        pre {
          color: var(--ms-text) !important;
          background: #17141B !important;
          border-color: var(--ms-border) !important;
        }
        .ms-canvas {
          position: relative;
          min-height: 520px;
          border: 1px solid var(--ms-line);
          border-radius: var(--ms-radius-lg);
          overflow: hidden;
          background:
            radial-gradient(circle at 44% 28%, rgba(255, 98, 180, 0.20), transparent 18rem),
            radial-gradient(circle at 70% 30%, rgba(140, 218, 223, 0.12), transparent 18rem),
            linear-gradient(180deg, rgba(21, 18, 27, 0.82), rgba(9, 8, 12, 0.94));
          box-shadow: inset 0 0 90px rgba(0,0,0,0.38), 0 18px 50px rgba(0,0,0,0.28);
        }
        .ms-legend {
          position: absolute;
          top: 1.1rem;
          left: 1.15rem;
          display: grid;
          gap: 0.45rem;
          z-index: 3;
        }
        .ms-legend span {
          display: flex;
          align-items: center;
          gap: 0.5rem;
          color: var(--ms-ink);
          font-weight: 700;
        }
        .ms-legend b {
          display: inline-block;
          width: 0.75rem;
          height: 0.75rem;
          border-radius: 50%;
        }
        .ms-legend .cyan { background: var(--ms-teal); }
        .ms-legend .amber { background: var(--ms-amber); }
        .ms-legend .silver { background: #aeb8c1; }
        .ms-protein-shape {
          position: absolute;
          left: 21%;
          right: 10%;
          top: 8%;
          height: 50%;
          filter: drop-shadow(0 28px 38px rgba(0,0,0,0.38));
        }
        .ms-protein-shape i {
          position: absolute;
          display: block;
          border-radius: 50%;
          background: radial-gradient(circle at 28% 25%, #a9f6ff, var(--ms-teal) 45%, #176b87);
          opacity: 0.94;
        }
        .ms-protein-shape i:nth-child(1) { width: 16rem; height: 13rem; left: 0; top: 2.2rem; }
        .ms-protein-shape i:nth-child(2) { width: 13rem; height: 12rem; left: 11rem; top: 0.5rem; background: radial-gradient(circle at 30% 25%, #c9bcff, #6967d8 50%, #343681); }
        .ms-protein-shape i:nth-child(3) { width: 12rem; height: 11rem; left: 5.5rem; top: 11rem; background: radial-gradient(circle at 30% 25%, #8bd8ff, #357de1 55%, #174b92); }
        .ms-protein-shape i:nth-child(4) { width: 13rem; height: 10.5rem; left: 21rem; top: 11.4rem; background: radial-gradient(circle at 30% 25%, #c5fff5, #61c8c0 52%, #286c69); }
        .ms-protein-shape i:nth-child(5) { width: 8rem; height: 8rem; left: 27rem; top: 4rem; background: radial-gradient(circle at 30% 25%, #ff9fca, var(--ms-accent) 55%, #9b1f61); }
        .ms-protein-shape i:nth-child(6) { width: 7rem; height: 7rem; left: 14.5rem; top: 8.5rem; background: radial-gradient(circle at 30% 25%, #ffffff, #76e4ef 52%, #2d7e95); }
        .ms-surface-grid {
          position: absolute;
          left: -4%;
          right: -4%;
          bottom: -7%;
          height: 34%;
          background:
            radial-gradient(circle, rgba(210, 222, 231, 0.36) 0 0.48rem, transparent 0.52rem) 0 0 / 2.1rem 1.35rem,
            radial-gradient(circle, rgba(130, 148, 160, 0.38) 0 0.46rem, transparent 0.52rem) 1.05rem 0.68rem / 2.1rem 1.35rem,
            linear-gradient(180deg, rgba(255,255,255,0.12), rgba(0,0,0,0.24));
          transform: perspective(580px) rotateX(58deg);
          transform-origin: bottom;
        }
        .ms-anchor-callout {
          position: absolute;
          bottom: 34%;
          border: 1px solid rgba(255, 185, 76, 0.55);
          background: rgba(5, 16, 28, 0.82);
          color: var(--ms-amber);
          border-radius: var(--ms-radius-md);
          padding: 0.45rem 0.65rem;
          font-weight: 800;
          z-index: 4;
        }
        .ms-anchor-callout.left { left: 9%; }
        .ms-anchor-callout.right { right: 8%; color: var(--ms-teal); border-color: rgba(140, 218, 223, 0.45); }
        .ms-axis {
          position: absolute;
          left: 2rem;
          bottom: 2rem;
          color: var(--ms-ink);
          font-weight: 900;
        }
        .ms-spec-grid {
          display: grid;
          grid-template-columns: repeat(4, minmax(0, 1fr));
          gap: 0.8rem;
          margin: 1rem 0;
        }
        .ms-spec {
          border: 1px solid var(--ms-line);
          border-radius: var(--ms-radius-lg);
          background: rgba(33, 28, 39, 0.82);
          padding: 0.95rem;
        }
        .ms-spec span {
          display: block;
          color: var(--ms-muted);
          font-size: 0.82rem;
          margin-bottom: 0.32rem;
        }
        .ms-spec strong {
          color: var(--ms-ink);
          font-size: 1rem;
        }
        .ms-command-bar {
          border: 1px solid var(--ms-line);
          border-radius: var(--ms-radius-lg);
          background: rgba(18, 16, 22, 0.82);
          padding: 0.6rem;
        }
        .ms-empty-preview {
          min-height: 430px;
          border: 1px dashed rgba(255, 98, 180, 0.28);
          border-radius: var(--ms-radius-lg);
          background:
            radial-gradient(circle at center, rgba(255, 92, 173, 0.08), transparent 18rem),
            rgba(18, 16, 22, 0.70);
          display: grid;
          place-content: center;
          text-align: center;
          padding: 2rem;
        }
        .ms-empty-preview strong {
          color: var(--ms-ink);
          font-size: 1.4rem;
        }
        .ms-empty-preview span {
          color: var(--ms-muted);
          max-width: 32rem;
          margin-top: 0.45rem;
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
          color: var(--ms-text);
          background: rgba(255, 98, 180, 0.12);
          padding: 0.28rem 0.68rem;
          font-size: 0.82rem;
          font-weight: 650;
        }
        .ms-section {
          border: 1px solid var(--ms-line);
          border-radius: var(--ms-radius-lg);
          background: linear-gradient(180deg, rgba(33, 28, 39, 0.90), rgba(18, 16, 22, 0.88));
          padding: 1rem;
          margin-bottom: 1rem;
          box-shadow: var(--ms-shadow-panel);
        }
        .ms-log-row {
          border-left: 3px solid var(--ms-accent);
          padding: 0.35rem 0 0.35rem 0.75rem;
          margin: 0.35rem 0;
          background: rgba(33, 28, 39, 0.72);
          border-radius: 0 var(--ms-radius-md) var(--ms-radius-md) 0;
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
          border-radius: var(--ms-radius-md);
          border: 1px solid var(--ms-border);
          background: rgba(27, 23, 32, 0.92);
          color: var(--ms-text);
          min-height: 2.55rem;
          font-weight: 780;
          box-shadow: 0 5px 16px rgba(0,0,0,0.12);
          transition: background-color 160ms ease, border-color 160ms ease, transform 160ms ease;
        }
        .stButton > button:hover,
        .stDownloadButton > button:hover {
          border-color: var(--ms-primary);
          background: #2B2432;
          color: var(--ms-text);
          transform: translateY(-1px);
        }
        .stButton > button:focus,
        .stDownloadButton > button:focus,
        .stLinkButton > a:focus {
          box-shadow: 0 0 0 0.18rem var(--ms-focus) !important;
          outline: none !important;
        }
        .stButton > button:disabled,
        .stDownloadButton > button:disabled {
          border-color: rgba(83, 70, 88, 0.55) !important;
          background: rgba(27, 23, 32, 0.46) !important;
          color: var(--ms-text-muted) !important;
          cursor: not-allowed !important;
        }
        .stLinkButton > a,
        div[data-testid="stPopover"] button,
        .stDownloadButton > button {
          border-radius: var(--ms-radius-md) !important;
          border: 1px solid var(--ms-border) !important;
          background: rgba(27, 23, 32, 0.86) !important;
          color: var(--ms-text) !important;
          font-weight: 760 !important;
          min-height: 2.55rem;
        }
        .stLinkButton > a p,
        div[data-testid="stPopover"] button p,
        .stDownloadButton > button p {
          color: var(--ms-text) !important;
          overflow-wrap: normal !important;
          word-break: normal !important;
          white-space: normal !important;
        }
        .stButton > button[kind="primary"] {
          background: linear-gradient(135deg, var(--ms-accent), #FF7BBC) !important;
          color: #180712 !important;
          border-color: rgba(255, 98, 180, 0.74) !important;
          font-weight: 880;
          box-shadow: 0 12px 28px rgba(255, 98, 180, 0.22);
        }
        .stButton > button[kind="primary"]:hover {
          background: linear-gradient(135deg, #FF69B1, #FF90C9) !important;
          color: #160711 !important;
        }
        div[data-testid="stAlert"] {
          border-radius: var(--ms-radius-md);
          border: 1px solid var(--ms-border);
          background: rgba(33, 28, 39, 0.94);
          color: var(--ms-text);
        }
        div[data-testid="stAlert"] p,
        div[data-testid="stAlert"] li {
          color: var(--ms-text) !important;
        }
        div[data-testid="stDataFrame"],
        div[data-testid="stDataEditor"] {
          border: 1px solid var(--ms-border);
          border-radius: var(--ms-radius-lg);
          overflow: hidden;
          background: rgba(24, 21, 29, 0.86);
          box-shadow: var(--ms-shadow-panel);
        }
        div[data-testid="stDataFrame"] *,
        div[data-testid="stDataEditor"] * {
          color: inherit;
        }
        div[data-testid="stStatusWidget"],
        div[data-testid="stStatus"] {
          color: var(--ms-text) !important;
        }
        @media (max-width: 900px) {
          .ms-header {
            grid-template-columns: 1fr;
          }
          .ms-topbar {
            grid-template-columns: 1fr;
          }
          .ms-hero-logo {
            display: none;
          }
          .ms-spec-grid {
            grid-template-columns: 1fr 1fr;
          }
          .ms-canvas {
            min-height: 420px;
          }
        }
        @media (max-width: 1100px) {
          .block-container {
            padding-left: 1.25rem;
            padding-right: 1.25rem;
          }
          section[data-testid="stMain"] div[data-testid="stHorizontalBlock"] {
            display: block !important;
          }
          section[data-testid="stMain"] div[data-testid="stHorizontalBlock"] > div[data-testid="stColumn"] {
            flex: 1 1 100% !important;
            min-width: 100% !important;
            width: 100% !important;
            margin-bottom: 0.75rem;
          }
          .stLinkButton > a,
          div[data-testid="stPopover"] button,
          .stDownloadButton > button,
          .stButton > button {
            min-height: 2.75rem;
          }
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
