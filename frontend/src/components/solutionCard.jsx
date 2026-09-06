import React, { useEffect, useRef, useState } from 'react';
import './solutionCard.css';

const solutionText = {
  title: "Revolutionize Biodiversity Discovery with AI",
  subtitle: "Transform environmental DNA samples into actionable biodiversity insights with cutting-edge AI technology.",
  bullets: [
    "Automated species identification from raw sequencing data",
    "Real-time quality control and contamination detection",
    "Interactive visualization of taxonomic abundance and diversity"
  ],
  secondaryLabel: "Learn More"
};

export default function SolutionCard({ onNavigate }) {
  const [isVisible, setIsVisible] = useState(false);
  const cardRef = useRef(null);

  useEffect(() => {
    const observer = new IntersectionObserver(
      ([entry]) => {
        if (entry.isIntersecting) {
          setIsVisible(true);
          observer.disconnect();
        }
      },
      { threshold: 0.1 }
    );

    if (cardRef.current) {
      observer.observe(cardRef.current);
    }

    return () => observer.disconnect();
  }, []);

  return (
    <section className="solution-section">
      <div
        ref={cardRef}
        className={`solution-card ${isVisible ? 'visible' : ''}`}
      >
        <div className="solution-content">
          {/* Top image */}
          <div className="solution-image-container animate-item">
            <div className="solution-image-wrapper">
              <img
                src="/assets/images/cards/Sol.jpg"
                alt="AI Biodiversity Discovery Solution"
                loading="lazy"
                className="solution-image"
              />
            </div>
          </div>

          {/* Text content */}
          <div className="solution-text">
            <h2 id="solution-title" className="solution-title animate-item">
              {solutionText.title}
            </h2>
            <p className="solution-subtitle animate-item">
              {solutionText.subtitle}
            </p>

            <ul className="solution-bullets animate-item">
              {solutionText.bullets.map((bullet, i) => (
                <li key={i} className="solution-bullet-item">
                  <span className="bullet-dot" aria-hidden="true"></span>
                  <span>{bullet}</span>
                </li>
              ))}
            </ul>

            <div className="solution-actions animate-item">
              <button 
                className="btn-secondary"
                onClick={() => onNavigate('about')}
              >
                {solutionText.secondaryLabel}
              </button>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
}