import React, { useEffect, useRef } from 'react';
import Typewriter from 'typewriter-effect';
import SolutionCard from '../components/solutionCard';
import WhyChoose from '../components/WhyChoose';
import FeaturesSection from '../components/FeaturesSection';
import Footer from '../components/Footer';

import './home.css';

export default function Home({ onNavigate }) {
  const heroRef = useRef(null);

  useEffect(() => {
    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            entry.target.classList.add('fade-in');
          }
        });
      },
      { threshold: 0.1 }
    );

    const elements = heroRef.current?.querySelectorAll('.animate-on-scroll');
    elements?.forEach((el) => observer.observe(el));

    return () => observer.disconnect();
  }, []);

  return (
    <main className="home-page">
      <section className="hero" ref={heroRef}>
        <div className="hero-content">
          <div className="hero-left animate-on-scroll">
            <h1 className="hero-title">
              AI-POWERED eDNA ANALYSIS PLATFORM
            </h1>
            <div className="hero-subtitle">
              <span className="subtitle-prefix">Transforming </span>
              <span className="typing-text">
                <Typewriter
                  options={{
                    strings: [
                      'upload',
                      'download',
                      'research',
                      'analyze',
                      'raw data',
                      'sequencing',
                      'reads',
                      'diversity insights at speed'
                    ],
                    autoStart: true,
                    loop: true,
                    delay: 75,
                    deleteSpeed: 50,
                    pauseFor: 2000,
                  }}
                />
              </span>
              <span className="subtitle-suffix"> into actionable insights</span>
            </div>
            <button 
              className="btn-primary hero-cta" 
              onClick={() => onNavigate('run')}
            >
              Start Analysis
            </button>
          </div>

          <div className="hero-right animate-on-scroll">
            <div className="video-container">
              <video
                className="hero-video"
                src="/assets/hero-video.mp4"
                muted
                loop
                autoPlay
                playsInline
                poster="/assets/hero-fallback.png"
              >
                Your browser does not support the video tag.
              </video>
            </div>
          </div>
        </div>
      </section>

      <SolutionCard onNavigate={onNavigate} />

      <WhyChoose />

      <FeaturesSection onNavigate={onNavigate} />

      <section className="features-preview animate-on-scroll">
        <h2>Accelerate Your Research</h2>
        <p>Process environmental DNA samples with cutting-edge AI technology</p>
      </section>

      <Footer />
    </main>
  );
}